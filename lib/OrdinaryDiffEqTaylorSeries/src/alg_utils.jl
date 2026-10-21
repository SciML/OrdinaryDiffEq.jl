get_value(::Val{T}) where {T} = T

alg_order(::ExplicitTaylor2) = 2
alg_stability_size(alg::ExplicitTaylor2) = 1

alg_order(::ExplicitTaylor{P}) where {P} = P
alg_stability_size(alg::ExplicitTaylor) = 1

# These build each step from the Taylor jet rather than reusing the previous step's
# last derivative, so they are not FSAL. Without this they inherit the `true` default,
# which makes `ode_determine_initdt` write an extra `f` evaluation into `fsallast`.
isfsal(::ExplicitTaylor) = false
isfsal(::ExplicitTaylorAdaptiveOrder) = false

alg_order(alg::ExplicitTaylorAdaptiveOrder) = get_value(alg.min_order)
get_current_adaptive_order(::ExplicitTaylorAdaptiveOrder, cache) = cache.current_order[]
get_current_alg_order(::ExplicitTaylorAdaptiveOrder, cache) = cache.current_order[]

# A step at `current_order` is jetted one order higher, so `current_order` has to
# stay below `max_order`. With an empty window the order loop never runs and the
# step is accepted on a stale error estimate.
function order_window(alg::ExplicitTaylorAdaptiveOrder)
    min_order_value = get_value(alg.min_order)
    max_order_value = get_value(alg.max_order)
    if max_order_value <= min_order_value
        throw(
            ArgumentError(
                "ExplicitTaylorAdaptiveOrder requires max_order > min_order, got " *
                    "min_order = $min_order_value, max_order = $max_order_value"
            )
        )
    end
    return min_order_value, max_order_value
end

JET_CACHE = IdDict()
# f -> (coeffs, params, polynomial, d_polynomial)
POLYNOMIAL_CACHE = IdDict()

# Opt-in joint evaluators used by the implicit solver.
POLYNOMIAL_BUNDLE_CACHE = IdDict()

# Symbolics' nested in-place output writes each buffer back to its container.
# Keep heterogeneous buffers concretely typed, while retaining their storage.
struct PolynomialOutputs{T}
    arrays::T
end
@inline Base.getindex(out::PolynomialOutputs, i::Int) = out.arrays[i]
@inline function Base.setindex!(out::PolynomialOutputs, value, i::Int)
    value === out.arrays[i] || copyto!(out.arrays[i], value)
    return value
end

# A joint kernel evaluating several weighted Taylor polynomials (and optionally
# the Jacobian of the first one) from a single jet at `(u, t)`. The outputs
# come first, then `u, t, h`.
struct PolynomialBundle{F}
    f::F
end
@inline function (bundle::PolynomialBundle)(a, u, t, h)
    bundle.f(PolynomialOutputs((a,)), u, t, h)
    return nothing
end
@inline function (bundle::PolynomialBundle)(a, b, u, t, h)
    bundle.f(PolynomialOutputs((a, b)), u, t, h)
    return nothing
end
@inline function (bundle::PolynomialBundle)(a, b, c, u, t, h)
    bundle.f(PolynomialOutputs((a, b, c)), u, t, h)
    return nothing
end
@inline function (bundle::PolynomialBundle)(a, b, c, d, u, t, h)
    bundle.f(PolynomialOutputs((a, b, c, d)), u, t, h)
    return nothing
end

# Helper functions to unwrap Symbolics.Num to concrete values (needed for Symbolics v7+)
# After build_function evaluation, results may be wrapped in Num and need explicit unwrapping

# Convert symbolic value to concrete numeric type
# Uses Symbolics.value for proper evaluation to concrete numbers
@inline function _sym_to_numeric(x::Symbolics.Num, ::Type{T}) where {T}
    v = SymbolicUtils.unwrap(x)
    return v isa Number ? convert(T, v) : convert(T, Symbolics.value(x))
end
@inline function _sym_to_numeric(x::SymbolicUtils.BasicSymbolic, ::Type{T}) where {T}
    # BasicSymbolic that should be a concrete number - evaluate it
    return convert(T, Symbolics.value(Symbolics.Num(x)))
end
@inline _sym_to_numeric(x::Number, ::Type{T}) where {T} = convert(T, x)
@inline _sym_to_numeric(x, ::Type{T}) where {T} = convert(T, x)  # fallback

# Unwrap Symbolics values in TaylorScalar coefficients
@inline function sym_unwrap_taylor(ts::TaylorScalar, ::Type{T}) where {T}
    coeffs = TaylorDiff.flatten(ts)
    unwrapped = map(c -> _sym_to_numeric(c, T), coeffs)
    return TaylorScalar(unwrapped)  # Pass tuple directly
end

# For arrays of TaylorScalars
@inline function sym_unwrap_taylor(arr::AbstractArray{<:TaylorScalar}, ::Type{T}) where {T}
    return map(ts -> sym_unwrap_taylor(ts, T), arr)
end

# Fallback for other types (e.g., already-concrete values)
@inline sym_unwrap_taylor(x, ::Type{T}) where {T} = x

# A jet of order `P` does not always fill a Taylor buffer of the same order:
# `ExplicitTaylorAdaptiveOrder` sizes one buffer at `max_order` and runs every
# order through it. Assigning a `TaylorScalar{T, P}` where a `TaylorScalar{T, Q}`
# is expected goes through TaylorDiff's number-promotion constructor, which keeps
# the constant term and silently drops every derivative, so widen through the
# order-changing constructor instead. It zero-pads coefficients above `P`, which
# is exactly the polynomial the order-`P` jet represents.
@inline set_taylor_order(t::TaylorScalar, ::Val{Q}) where {Q} = TaylorScalar{Q}(t)
# `ExplicitTaylor` always jets straight into a buffer of its own order; rebuilding the
# `TaylorScalar` there is a copy that the compiler does not always elide.
@inline set_taylor_order(t::TaylorScalar{T, Q}, ::Val{Q}) where {T, Q} = t

@inline function taylor_from_coefficients(coeffs, i, ::Val{P}, buffer_order) where {P}
    offset = (i - 1) * (P + 1)
    taylor = TaylorScalar(ntuple(j -> coeffs[offset + j], Val(P + 1)))
    return set_taylor_order(taylor, buffer_order)
end

function build_jet(
        f::ODEFunction{iip}, p, order, length = nothing, buffer_order = order
    ) where {iip}
    # Unwrap FunctionWrappers since TaylorDiff/Symbolics types don't match wrapper signatures
    f = unwrapped_f(f)
    return build_jet(f, Val{iip}(), p, order, length, buffer_order)
end

# cache format: Dict{typeof(f), Vector{Tuple{order, p, buffer_order, jet}}}
function build_jet(
        f, ::Val{iip}, p, order::Val{P}, length = nothing,
        buffer_order::Val{Q} = order
    ) where {P, Q, iip}
    if haskey(JET_CACHE, f)
        list = JET_CACHE[f]
        index = findfirst(x -> x[1] == order && x[2] == p && x[3] == buffer_order, list)
        index !== nothing && return list[index][4]
    end
    @variables t0::Real
    u0 = isnothing(length) ? Symbolics.variable(:u0) : Symbolics.variables(:u0, 1:length)
    if iip
        @assert length isa Integer
        f0 = similar(u0)
        f(f0, u0, p, t0)
    else
        f0 = f(u0, p, t0)
    end
    u = TaylorDiff.make_seed(u0, f0, Val(1))
    for index in 2:P
        t = TaylorScalar{index - 1}(t0, one(t0))
        if iip
            fu = similar(u)
            f(fu, u, p, t)
        else
            fu = f(u, p, t)
        end
        d = TaylorDiff.get_coefficient(fu, index - 1) / index
        u = TaylorDiff.append_coefficient(u, d)
    end

    # Flatten TaylorScalar coefficients for build_function (it doesn't handle custom structs)
    # Then wrap result back into TaylorScalar
    if u isa TaylorScalar
        # Scalar case: build function for array of coefficients
        coeffs = collect(TaylorDiff.flatten(u))
        scalar_jet_coeffs = build_function(
            coeffs, u0, t0; expression = Val(false), cse = true
        )
        # Wrap to return TaylorScalar
        jet = (u0_val, t0_val) -> set_taylor_order(
            TaylorScalar(Tuple(scalar_jet_coeffs[1](u0_val, t0_val))), buffer_order
        )
    elseif u isa AbstractArray && eltype(u) <: TaylorScalar
        n = Base.length(u)
        coeffs = [TaylorDiff.flatten(u[i])[j] for i in 1:n for j in 1:(P + 1)]
        array_jet_coeffs = build_function(
            coeffs, u0, t0; expression = Val(false), cse = true
        )
        jet = (
            (u0_val, t0_val) -> begin
                coeffs_out = array_jet_coeffs[1](u0_val, t0_val)
                return [taylor_from_coefficients(coeffs_out, i, order, buffer_order) for i in 1:n]
            end,
            (out, coeffs_out, u0_val, t0_val) -> begin
                array_jet_coeffs[2](coeffs_out, u0_val, t0_val)
                for i in 1:n
                    out[i] = taylor_from_coefficients(coeffs_out, i, order, buffer_order)
                end
                return out
            end,
        )
    else
        # Fallback (shouldn't happen normally)
        jet = build_function(u, u0, t0; expression = Val(false), cse = true)
    end

    if !haskey(JET_CACHE, f)
        JET_CACHE[f] = []
    end
    push!(JET_CACHE[f], (order, p, buffer_order, jet))
    return jet
end

function build_polynomial(f::ODEFunction{iip}, p, coeffs::NTuple{P1, Float64}, length = nothing) where {P1, iip}
    f = unwrapped_f(f)
    return build_polynomial(f, Val{iip}(), p, coeffs, length)
end

function build_polynomial(f, ::Val{iip}, p, coeffs::NTuple{P1, Float64}, length = nothing) where {P1, iip}
    if haskey(POLYNOMIAL_CACHE, f)
        list = POLYNOMIAL_CACHE[f]
        index = findfirst(x -> x[1] == coeffs && x[2] == p, list)
        index !== nothing && return list[index][3], list[index][4]
    end
    expressions, u0, t0, dt, _ = polynomial_expressions(f, Val(iip), p, (coeffs,), length)
    ut = only(expressions)
    polynomial = build_function(ut, u0, t0, dt; expression = Val(false), cse = true)
    jacobian = Symbolics.jacobian(ut, u0)
    d_polynomial = build_function(jacobian, u0, t0, dt; expression = Val(false), cse = true)

    if !haskey(POLYNOMIAL_CACHE, f)
        POLYNOMIAL_CACHE[f] = []
    end
    push!(POLYNOMIAL_CACHE[f], (coeffs, p, polynomial, d_polynomial))
    return polynomial, d_polynomial
end

# Build the Taylor coefficients once, including total time derivatives, then
# form any number of weighted polynomials at the same state and time.
function polynomial_expressions(f, ::Val{iip}, p, coefficient_sets, length) where {iip}
    P = maximum(Base.length, coefficient_sets) - 1
    @variables t0::Real dt::Real
    u0 = isnothing(length) ? Symbolics.variable(:u0) : Symbolics.variables(:u0, 1:length)
    if iip
        f0 = similar(u0)
        f(f0, u0, p, t0)
    else
        f0 = f(u0, p, t0)
    end
    u = TaylorDiff.make_seed(u0, f0, Val(1))
    for index in 2:P
        t = TaylorScalar{index - 1}(t0, one(t0))
        if iip
            fu = similar(u)
            f(fu, u, p, t)
        else
            fu = f(u, p, t)
        end
        d = TaylorDiff.get_coefficient(fu, index - 1) / index
        u = TaylorDiff.append_coefficient(u, d)
    end
    expressions = map(coefficient_sets) do coeffs
        padded = ntuple(k -> k <= Base.length(coeffs) ? coeffs[k] : 0.0, P + 1)
        eval_taylor_polynomial(u, padded, dt)
    end
    return expressions, u0, t0, dt, u
end

# Jacobians ∂c_k/∂u of the Taylor coefficients c_k = y^(k)/k!, k = 1:Q, of
# a jet `u`, stacked into one n×n×Q array so that the generated kernel
# writes them in place as a single output. They do not depend on the step
# size, so a Newton matrix I + Σ_k w_k h^k ∂c_k/∂u can be assembled for any
# h from one evaluation.
function coefficient_jacobians_of(u, u0, Q)
    jacobians = [
        Symbolics.jacobian([TaylorDiff.flatten(u[i])[k + 1] for i in eachindex(u)], u0)
            for k in 1:Q
    ]
    return cat(jacobians...; dims = 3)
end

function build_polynomial_bundle(
        f::ODEFunction{iip}, p, coefficient_sets::Tuple, length;
        jacobian = false, coefficient_jacobians = 0
    ) where {iip}
    return build_polynomial_bundle(
        unwrapped_f(f), Val(iip), p, coefficient_sets,
        length; jacobian, coefficient_jacobians
    )
end

# Build one in-place kernel `bundle(out_1, ..., out_n[, J], u, t, h)` that
# evaluates every weighted polynomial in `coefficient_sets` at the same state
# and time from a single Taylor jet. With `jacobian = true` the Jacobian of
# the first polynomial with respect to `u` is appended as the last output;
# with `coefficient_jacobians = Q` the Jacobians of the Taylor coefficients
# c_1, …, c_Q are appended instead as one n×n×Q output. Joint code generation
# shares common expressions across the polynomial and Jacobian outputs.
function build_polynomial_bundle(
        f, ::Val{iip}, p, coefficient_sets::Tuple, length;
        jacobian = false, coefficient_jacobians = 0
    ) where {iip}
    kernels = get!(Dict{Tuple, Vector{Any}}, POLYNOMIAL_BUNDLE_CACHE, f)
    key = (coefficient_sets, length, jacobian, coefficient_jacobians, iip)
    entries = get!(Vector{Any}, kernels, key)
    index = findfirst(entry -> entry[1] == p, entries)
    index === nothing || return entries[index][2]
    expressions, u0, t0, dt, jet = polynomial_expressions(
        f, Val(iip), p, coefficient_sets, length
    )
    outputs = if jacobian
        [expressions..., Symbolics.jacobian(first(expressions), u0)]
    elseif coefficient_jacobians > 0
        [expressions..., coefficient_jacobians_of(jet, u0, coefficient_jacobians)]
    else
        collect(expressions)
    end
    # Nested outputs write directly into preallocated buffers. A single CSE pass
    # shares expressions across all outputs without copying a combined matrix.
    _, joint! = build_function(
        outputs, u0, t0, dt;
        expression = Val(false), cse = true, skipzeros = false
    )
    bundle = PolynomialBundle(joint!)
    push!(entries, (p, bundle))
    return bundle
end

# evaluate using Qin Jiushao's algorithm
@generated function evaluate_polynomial(t::TaylorScalar{T, P}, z) where {T, P}
    ex = :(v[$(P + 1)])
    for i in P:-1:1
        ex = :(v[$i] + z * $ex)
    end
    return :($(Expr(:meta, :inline)); v = TaylorDiff.flatten(t); $ex)
end

# Evaluate polynomial for scalar TaylorScalar (returns scalar)
@inline eval_taylor_polynomial(utaylor::TaylorScalar, dt) = evaluate_polynomial(utaylor, dt)
# Evaluate polynomial for array of TaylorScalars (returns array)
@inline eval_taylor_polynomial(utaylor::AbstractArray, dt) = map(x -> evaluate_polynomial(x, dt), utaylor)

# Evaluate a weighted Taylor polynomial for scalar or array states.
@inline eval_taylor_polynomial(u::TaylorScalar, coeffs, dt) = evalpoly(dt, map(*, coeffs, TaylorDiff.flatten(u)))
@inline eval_taylor_polynomial(us::AbstractArray, coeffs, dt) = map(x -> eval_taylor_polynomial(x, coeffs, dt), us)
