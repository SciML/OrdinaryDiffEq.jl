import OrdinaryDiffEqTaylorSeries: build_jet, build_polynomial

@cache mutable struct ImplicitTaylorCache{
        T, uType, tType, jacobianType, uIntermediateType, rateType, uNoUnitsType, F1, Tol, StepLimiter,
    } <:
    OrdinaryDiffEqMutableCache
    μ::T
    t::tType
    polynomial::FunctionWrapper{Nothing, Tuple{uType, uType, tType, tType}}
    d_polynomial::FunctionWrapper{Nothing, Tuple{jacobianType, uType, tType, tType}}
    polynomial_explicit::FunctionWrapper{Nothing, Tuple{uType, uType, tType, tType}}
    u::uType
    uprev::uType
    uprev2::uType
    J::jacobianType
    utilde::uType
    uintermediate::uIntermediateType
    rhs::uType
    tmp::uIntermediateType
    atmp::uNoUnitsType
    fsalfirst::rateType
    linsolve::F1
    κ::Tol
    ηold::Tol
    status::NLStatus
    iter::Int
    step_limiter!::StepLimiter
end

function alg_cache(
        alg::ImplicitTaylor{P, Q}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {P, Q, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if is_mu_taylor(alg)
        coeffs = ntuple(_ -> 1.0, alg.order)
        polynomial_all, d_polynomial_all = build_polynomial(f, p, coeffs, length(u))
        # get the iip version
        _, polynomial = polynomial_all
        _, d_polynomial = d_polynomial_all
        polynomial_explicit = polynomial
    else
        polynomial_p, polynomial_q = normalized_pade(P, Q)
        tuple_p = Base.tail(tuple(map(Float64, polynomial_p)...))
        tuple_q = Base.tail(tuple(map(Float64, polynomial_q)...))
        polynomial_p_all, _ = build_polynomial(f, p, tuple_p, length(u))
        polynomial_q_all, d_polynomial_q_all = build_polynomial(f, p, tuple_q, length(u))
        # get the iip version
        _, polynomial = polynomial_q_all
        _, d_polynomial = d_polynomial_q_all
        _, polynomial_explicit = polynomial_p_all
    end
    fsalfirst = zero(rate_prototype)
    atmp = similar(u, promote_type(uEltypeNoUnits, typeof(alg.μ)))
    recursivefill!(atmp, false)
    utilde = zero(u)
    uToltype = constvalue(uBottomEltypeNoUnits)
    promoted_eltype = promote_type(eltype(u), typeof(alg.μ))

    tmp = similar(u, promoted_eltype)
    uintermediate = similar(u, promoted_eltype)
    rhs = similar(u)
    tmp_vec = _vec(tmp)
    linsolve_u0 = copy(tmp_vec)
    κ = alg.κ !== nothing ? convert(uToltype, alg.κ) : convert(uToltype, 1 // 100)
    J = ArrayInterface.zeromatrix(tmp_vec)
    linprob = LinearProblem(J, tmp_vec; u0 = linsolve_u0)
    linsolve = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        assumptions = LinearSolve.OperatorAssumptions(true)
    )

    ηold = one(uToltype)
    iter = 10000
    status = Convergence
    uType, tType, jacobianType = typeof(u), typeof(t), typeof(J)
    polynomial_wrapped = FunctionWrapper{Nothing, Tuple{uType, uType, tType, tType}}(polynomial)
    d_polynomial_wrapped = FunctionWrapper{Nothing, Tuple{jacobianType, uType, tType, tType}}(d_polynomial)
    polynomial_explicit_wrapped = FunctionWrapper{Nothing, Tuple{uType, uType, tType, tType}}(polynomial_explicit)

    return ImplicitTaylorCache(
        alg.μ, t, polynomial_wrapped, d_polynomial_wrapped, polynomial_explicit_wrapped, u, uprev, uprev2, J, utilde, uintermediate, rhs, tmp, atmp, fsalfirst, linsolve,
        κ, ηold, status, iter, alg.step_limiter!
    )
end

function get_fsalfirstlast(cache::ImplicitTaylorCache, u)
    return (cache.fsalfirst, cache.fsalfirst)
end

mutable struct ImplicitTaylorConstantCache{N} <: OrdinaryDiffEqConstantCache
    nlsolver::N
end

function alg_cache(
        alg::ImplicitTaylor, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false)
    )
    return ImplicitTaylorConstantCache(nlsolver)
end

@inline function SciMLBase.get_tmp_cache(
        integrator, alg::Union{ImplicitTaylor},
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.tmp, cache.atmp)
end
