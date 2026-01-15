@cache mutable struct IRKCConstantCache{uType, rateType, N} <: OrdinaryDiffEqConstantCache
    minm::Int
    zprev::uType
    nlsolver::N
    du₁::rateType
    du₂::rateType
end

@cache mutable struct IRKCCache{uType, rateType, uNoUnitsType, N, C <: IRKCConstantCache} <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    gprev::uType
    gprev2::uType
    fsalfirst::rateType
    f1ⱼ₋₁::rateType
    f1ⱼ₋₂::rateType
    f2ⱼ₋₁::rateType
    atmp::uNoUnitsType
    nlsolver::N
    du₁::rateType
    du₂::rateType
    constantcache::C
end

function get_fsalfirstlast(cache::IRKCCache, u)
    return (cache.fsalfirst, du_alias_or_new(cache.nlsolver, cache.fsalfirst))
end

function alg_cache(
        alg::IRKC, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1.0, 1.0
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    zprev = u
    du₁ = rate_prototype
    du₂ = rate_prototype
    return IRKCConstantCache(50, zprev, nlsolver, du₁, du₂)
end

function alg_cache(
        alg::IRKC, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1.0, 1.0
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )

    gprev = zero(u)
    gprev2 = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    zprev = zero(u)
    f1ⱼ₋₁ = zero(rate_prototype)
    f1ⱼ₋₂ = zero(rate_prototype)
    f2ⱼ₋₁ = zero(rate_prototype)
    du₁ = zero(rate_prototype)
    du₂ = zero(rate_prototype)
    constantcache = IRKCConstantCache(50, zprev, nlsolver, du₁, du₂)
    return IRKCCache(
        u, uprev, gprev, gprev2, fsalfirst, f1ⱼ₋₁, f1ⱼ₋₂, f2ⱼ₋₁, atmp, nlsolver, du₁,
        du₂, constantcache
    )
end
