import OrdinaryDiffEqTaylorSeries: build_jet, build_propagator

@cache mutable struct ImplicitTaylorCache{
        T, uType, propagatorType, d_propagatorType, jacobianType, rateType, uNoUnitsType, F1, Tol, StepLimiter,
    } <:
    OrdinaryDiffEqMutableCache
    μ::T
    propagator::propagatorType
    d_propagator::d_propagatorType
    u::uType
    uprev::uType
    uprev2::uType
    ut::uType
    J::jacobianType
    utilde::uType
    tmp::uType
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
        alg::ImplicitTaylor, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    propagator_all, d_propagator_all = build_propagator(f, p, alg.order, length(u))
    # get the iip version
    _, propagator = propagator_all
    _, d_propagator = d_propagator_all
    fsalfirst = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    utilde = zero(u)
    tmp = zero(u)
    uToltype = constvalue(uBottomEltypeNoUnits)
    ut = zero(u)

    κ = alg.κ !== nothing ? convert(uToltype, alg.κ) : convert(uToltype, 1 // 100)
    J = ArrayInterface.zeromatrix(_vec(u))
    linprob = LinearProblem(J, _vec(tmp); u0 = _vec(copy(u)))
    linsolve = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        assumptions = LinearSolve.OperatorAssumptions(true)
    )

    ηold = one(uToltype)
    iter = 10000
    status = Convergence

    return ImplicitTaylorCache(
        alg.μ, propagator, d_propagator, u, uprev, uprev2, ut, J, utilde, tmp, atmp, fsalfirst, linsolve,
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
