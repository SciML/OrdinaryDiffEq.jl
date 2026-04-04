get_value(::Val{T}) where {T} = T

alg_order(::ExplicitTaylor2) = 2
alg_stability_size(alg::ExplicitTaylor2) = 1

alg_order(::ExplicitTaylor{P}) where {P} = P
alg_stability_size(alg::ExplicitTaylor) = 1

alg_order(alg::ExplicitTaylorAdaptiveOrder) = get_value(alg.min_order)
get_current_adaptive_order(::ExplicitTaylorAdaptiveOrder, cache) = cache.current_order[]
get_current_alg_order(::ExplicitTaylorAdaptiveOrder, cache) = cache.current_order[]

# f -> (order, params, jet)
JET_CACHE = IdDict()
# f -> (coeffs, params, polynomial, d_polynomial)
POLYNOMIAL_CACHE = IdDict()

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

function build_polynomial(f::ODEFunction{iip}, p, coeffs::NTuple{P1, Float64}, length = nothing) where {P1, iip}
    f = unwrapped_f(f)
    return build_polynomial(f, Val{iip}(), p, coeffs, length)
end

function build_polynomial(f, ::Val{iip}, p, coeffs::NTuple{P1, Float64}, length = nothing) where {P1, iip}
    P = P1 - 1
    if haskey(POLYNOMIAL_CACHE, f)
        list = POLYNOMIAL_CACHE[f]
        index = findfirst(x -> x[1] == coeffs && x[2] == p, list)
        index !== nothing && return list[index][3:4]
    end
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
        d = get_coefficient(fu, index - 1) / index
        u = append_coefficient(u, d)
    end
    ut = eval_taylor_polynomial(u, coeffs, dt)
    polynomial = build_function(ut, u0, t0, dt; expression = Val(false), cse = true)
    jacobian = Symbolics.jacobian(ut, u0)
    d_polynomial = build_function(jacobian, u0, t0, dt; expression = Val(false), cse = true)

    if !haskey(POLYNOMIAL_CACHE, f)
        POLYNOMIAL_CACHE[f] = []
    end
    push!(POLYNOMIAL_CACHE[f], (coeffs, p, polynomial, d_polynomial))
    return polynomial, d_polynomial
end

# Evaluate polynomial for scalar TaylorScalar (returns scalar)
@inline eval_taylor_polynomial(u::TaylorScalar, coeffs, dt) = evalpoly(dt, map(*, coeffs, TaylorDiff.flatten(u)))
# Evaluate polynomial for array of TaylorScalars (returns array)
@inline eval_taylor_polynomial(us::AbstractArray, coeffs, dt) = map(x -> eval_taylor_polynomial(x, coeffs, dt), us)
