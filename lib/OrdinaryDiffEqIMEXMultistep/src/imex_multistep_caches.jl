# IMEX Multistep methods
abstract type IMEXMutableCache <: OrdinaryDiffEqMutableCache end
function get_fsalfirstlast(cache::IMEXMutableCache, u)
    return (cache.fsalfirst, du_alias_or_new(cache.nlsolver, cache.fsalfirst))
end

# CNAB2

@cache mutable struct CNAB2ConstantCache{rateType, N, uType, tType} <:
    OrdinaryDiffEqConstantCache
    k2::rateType
    nlsolver::N
    uprev3::uType
    tprev2::tType
end

@cache mutable struct CNAB2Cache{uType, rateType, N, tType} <: IMEXMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    fsalfirst::rateType
    k1::rateType
    k2::rateType
    du₁::rateType
    nlsolver::N
    uprev3::uType
    tprev2::tType
end

function alg_cache(
        alg::CNAB2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1 // 2, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )

    k2 = rate_prototype
    uprev3 = u
    tprev2 = t

    return CNAB2ConstantCache(k2, nlsolver, uprev3, tprev2)
end

function alg_cache(
        alg::CNAB2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1 // 2, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    du₁ = zero(rate_prototype)
    uprev3 = zero(u)
    tprev2 = t

    return CNAB2Cache(u, uprev, uprev2, fsalfirst, k1, k2, du₁, nlsolver, uprev3, tprev2)
end

# CNLF2

@cache mutable struct CNLF2ConstantCache{rateType, N, uType, tType} <:
    OrdinaryDiffEqConstantCache
    k2::rateType
    nlsolver::N
    uprev2::uType
    uprev3::uType
    tprev2::tType
end

@cache mutable struct CNLF2Cache{uType, rateType, N, tType} <: IMEXMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    fsalfirst::rateType
    k1::rateType
    k2::rateType
    du₁::rateType
    nlsolver::N
    uprev3::uType
    tprev2::tType
end

function alg_cache(
        alg::CNLF2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1 // 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )

    k2 = rate_prototype
    uprev2 = u
    uprev3 = u
    tprev2 = t

    return CNLF2ConstantCache(k2, nlsolver, uprev2, uprev3, tprev2)
end

function alg_cache(
        alg::CNLF2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1 // 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    du₁ = zero(rate_prototype)
    uprev2 = zero(u)
    uprev3 = zero(u)
    tprev2 = t

    return CNLF2Cache(u, uprev, uprev2, fsalfirst, k1, k2, du₁, nlsolver, uprev3, tprev2)
end
