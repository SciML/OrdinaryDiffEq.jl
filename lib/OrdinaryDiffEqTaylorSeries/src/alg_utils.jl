alg_order(::ExplicitTaylor2) = 2
alg_stability_size(alg::ExplicitTaylor2) = 1

alg_order(::ExplicitTaylor{P}) where {P} = P
alg_stability_size(alg::ExplicitTaylor) = 1

JET_CACHE = IdDict()

# Helper functions to unwrap Symbolics.Num to concrete values (needed for Symbolics v7+)
# After build_function evaluation, results may be wrapped in Num and need explicit unwrapping

# Convert symbolic value to concrete numeric type
# Uses Symbolics.value for proper evaluation to concrete numbers
@inline function _sym_to_numeric(x::Symbolics.Num, ::Type{T}) where {T}
    v = Symbolics.unwrap(x)
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

function build_jet(f::ODEFunction{iip}, p, order, length = nothing) where {iip}
    # Unwrap FunctionWrappers since TaylorDiff/Symbolics types don't match wrapper signatures
    f = unwrapped_f(f)
    return build_jet(f, Val{iip}(), p, order, length)
end

function build_jet(f, ::Val{iip}, p, order::Val{P}, length = nothing) where {P, iip}
    if haskey(JET_CACHE, f)
        list = JET_CACHE[f]
        index = findfirst(x -> x[1] == order && x[2] == p, list)
        index !== nothing && return list[index][3]
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
        d = get_coefficient(fu, index - 1) / index
        u = append_coefficient(u, d)
    end

    # Flatten TaylorScalar coefficients for build_function (it doesn't handle custom structs)
    # Then wrap result back into TaylorScalar
    if u isa TaylorScalar
        # Scalar case: build function for array of coefficients
        coeffs = collect(TaylorDiff.flatten(u))
        jet_coeffs = build_function(coeffs, u0, t0; expression = Val(false), cse = true)
        # Wrap to return TaylorScalar
        jet = (u0_val, t0_val) -> TaylorScalar(Tuple(jet_coeffs[1](u0_val, t0_val)))
    elseif u isa AbstractArray && eltype(u) <: TaylorScalar
        # Array case: build function for matrix of coefficients
        # Each row is the coefficients of one TaylorScalar
        n = Base.length(u)
        coeffs_matrix = [TaylorDiff.flatten(u[i])[j] for i in 1:n, j in 1:(P + 1)]
        jet_coeffs = build_function(coeffs_matrix, u0, t0; expression = Val(false), cse = true)
        # Wrap to return array of TaylorScalars
        jet = (
            (u0_val, t0_val) -> begin
                coeffs_out = jet_coeffs[1](u0_val, t0_val)
                return [TaylorScalar(Tuple(coeffs_out[i, :])) for i in 1:n]
            end,
            (out, u0_val, t0_val) -> begin
                coeffs_out = jet_coeffs[2](
                    similar(coeffs_matrix, eltype(u0_val)), u0_val, t0_val
                )
                for i in 1:n
                    out[i] = TaylorScalar(Tuple(coeffs_out[i, :]))
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
    push!(JET_CACHE[f], (order, p, jet))
    return jet
end

# evaluate using Qin Jiushao's algorithm
@generated function evaluate_polynomial(t::TaylorScalar{T, P}, z) where {T, P}
    ex = :(v[$(P + 1)])
    for i in P:-1:1
        ex = :(v[$i] + z * $ex)
    end
    return :($(Expr(:meta, :inline)); v = flatten(t); $ex)
end

# Evaluate polynomial for scalar TaylorScalar (returns scalar)
@inline eval_taylor_polynomial(utaylor::TaylorScalar, dt) = evaluate_polynomial(utaylor, dt)
# Evaluate polynomial for array of TaylorScalars (returns array)
@inline eval_taylor_polynomial(utaylor::AbstractArray, dt) = map(x -> evaluate_polynomial(x, dt), utaylor)
