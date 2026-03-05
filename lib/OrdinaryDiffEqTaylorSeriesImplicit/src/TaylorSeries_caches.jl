import OrdinaryDiffEqTaylorSeries: build_jet, build_propagator

abstract type ImplicitTaylorMutableCache <: OrdinaryDiffEqMutableCache end
abstract type ImplicitTaylorConstantCache <: OrdinaryDiffEqConstantCache end
function get_fsalfirstlast(cache::ImplicitTaylorMutableCache, u)
    return (cache.fsalfirst, cache.fsalfirst)
end

@cache mutable struct ImplicitTaylor1Cache{
        T, uType, propagatorType, d_propagatorType, jacobianType, rateType, uNoUnitsType, F1, Tol, StepLimiter,
    } <:
    ImplicitTaylorMutableCache
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
        alg::ImplicitTaylor1, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    propagator_all, d_propagator_all = build_propagator(f, p, Val(1), length(u))
    propagator, _ = propagator_all
    d_propagator, _ = d_propagator_all
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

    return ImplicitTaylor1Cache(
        alg.μ, propagator, d_propagator, u, uprev, uprev2, ut, J, utilde, tmp, atmp, fsalfirst, linsolve,
        κ, ηold, status, iter, alg.step_limiter!
    )
end

mutable struct ImplicitTaylor1ConstantCache{N} <: ImplicitTaylorConstantCache
    nlsolver::N
end

function alg_cache(
        alg::ImplicitTaylor1, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false)
    )
    return ImplicitTaylor1ConstantCache(nlsolver)
end

@inline function SciMLBase.get_tmp_cache(
        integrator, alg::Union{ImplicitTaylor1},
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.tmp, cache.atmp)
end
