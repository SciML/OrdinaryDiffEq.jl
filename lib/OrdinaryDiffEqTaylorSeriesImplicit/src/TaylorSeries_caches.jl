import OrdinaryDiffEqTaylorSeries: build_jet

abstract type ImplicitTaylorMutableCache <: OrdinaryDiffEqMutableCache end
abstract type ImplicitTaylorConstantCache <: OrdinaryDiffEqConstantCache end
function get_fsalfirstlast(cache::ImplicitTaylorMutableCache, u)
    (cache.fsalfirst, cache.fsalfirst)
end

@cache mutable struct ImplicitTaylor1Cache{
    T, uType, jetType, taylorType, rateType, uNoUnitsType, StepLimiter} <:
                      ImplicitTaylorMutableCache
    μ::T
    jet::jetType
    u::uType
    uprev::uType
    uprev2::uType
    utaylor::taylorType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    step_limiter!::StepLimiter
end

function alg_cache(alg::ImplicitTaylor1, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    jet_oop, jet_iip = build_jet(f, p, 1, length(u))
    utaylor = TaylorDiff.make_seed(u, zero(u), Val(1))
    # jet_wrapped = FunctionWrapper{typeof(utaylor), Tuple{typeof(u), typeof(t)}}(jet_oop)
    fsalfirst = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    utilde = zero(u)
    tmp = zero(u)

    ImplicitTaylor1Cache(
        alg.μ, jet_oop, u, uprev, uprev2, utaylor, utilde, tmp, atmp, fsalfirst,
        alg.step_limiter!)
end

mutable struct ImplicitTaylor1ConstantCache{N} <: ImplicitTaylorConstantCache
    nlsolver::N
end

function alg_cache(alg::ImplicitTaylor1, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1, 1
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false))
    ImplicitTaylor1ConstantCache(nlsolver)
end
