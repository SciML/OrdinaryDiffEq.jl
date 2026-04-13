abstract type StabilizedRKMutableCache <: OrdinaryDiffEqMutableCache end
get_fsalfirstlast(cache::StabilizedRKMutableCache, u) = (cache.fsalfirst, cache.k)

mutable struct ROCK2ConstantCache{T, T2, zType} <: OrdinaryDiffEqConstantCache
    ms::SVector{46, Int}
    fp1::SVector{46, T}
    fp2::SVector{46, T}
    recf::Vector{T2}
    zprev::zType
    mdeg::Int
    deg_index::Int
    start::Int
    min_stage::Int
    max_stage::Int
end
@cache struct ROCK2Cache{uType, rateType, uNoUnitsType, C <: ROCK2ConstantCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    k::rateType
    constantcache::C
end

function alg_cache(
        alg::ROCK2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    constantcache = ROCK2ConstantCache(
        constvalue(uBottomEltypeNoUnits),
        constvalue(tTypeNoUnits),
        u
    )
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return ROCK2Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::ROCK2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return ROCK2ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits), u)
end

mutable struct ROCK4ConstantCache{T, T2, T3, T4, zType} <: OrdinaryDiffEqConstantCache
    ms::SVector{50, Int}
    fpa::Vector{T}
    fpb::Vector{T2}
    fpbe::Vector{T3}
    recf::Vector{T4}
    zprev::zType
    mdeg::Int
    deg_index::Int
    start::Int
    min_stage::Int
    max_stage::Int
end

@cache struct ROCK4Cache{uType, rateType, uNoUnitsType, C <: ROCK4ConstantCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    uᵢ₋₃::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    k::rateType
    constantcache::C
end

function alg_cache(
        alg::ROCK4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    constantcache = ROCK4ConstantCache(
        constvalue(uBottomEltypeNoUnits),
        constvalue(tTypeNoUnits),
        u
    )
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    uᵢ₋₃ = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return ROCK4Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, uᵢ₋₃, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::ROCK4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return ROCK4ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits), u)
end

mutable struct RKCConstantCache{zType} <: OrdinaryDiffEqConstantCache
    #to match the types to call maxeig!
    zprev::zType
end
@cache struct RKCCache{uType, rateType, uNoUnitsType, C <: RKCConstantCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    gprev::uType
    gprev2::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    k::rateType
    constantcache::C
end

function alg_cache(
        alg::RKC, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    constantcache = RKCConstantCache(u)
    gprev = zero(u)
    gprev2 = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return RKCCache(u, uprev, gprev, gprev2, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::RKC, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RKCConstantCache(u)
end

mutable struct ESERK4ConstantCache{T, zType} <: OrdinaryDiffEqConstantCache
    ms::SVector{46, Int}
    Cᵤ::SVector{4, Int}
    Cₑ::SVector{4, Int}
    zprev::zType
    Bᵢ::Vector{T}
    mdeg::Int
    start::Int
    internal_deg::Int
end

@cache struct ESERK4Cache{uType, rateType, uNoUnitsType, C <: ESERK4ConstantCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Sᵢ::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    k::rateType
    constantcache::C
end

function alg_cache(
        alg::ESERK4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    constantcache = ESERK4ConstantCache(u)
    uᵢ = zero(u)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    Sᵢ = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return ESERK4Cache(u, uprev, uᵢ, uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::ESERK4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return ESERK4ConstantCache(u)
end

mutable struct ESERK5ConstantCache{T, zType} <: OrdinaryDiffEqConstantCache
    ms::SVector{49, Int}
    Cᵤ::SVector{5, Int}
    Cₑ::SVector{5, Int}
    zprev::zType
    Bᵢ::Vector{T}
    mdeg::Int
    start::Int
    internal_deg::Int
end

@cache struct ESERK5Cache{uType, rateType, uNoUnitsType, C <: ESERK5ConstantCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Sᵢ::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    k::rateType
    constantcache::C
end

function alg_cache(
        alg::ESERK5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    constantcache = ESERK5ConstantCache(u)
    uᵢ = zero(u)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    Sᵢ = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return ESERK5Cache(u, uprev, uᵢ, uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::ESERK5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return ESERK5ConstantCache(u)
end

mutable struct SERK2ConstantCache{T, zType} <: OrdinaryDiffEqConstantCache
    ms::SVector{11, Int}
    zprev::zType
    Bᵢ::Vector{T}
    mdeg::Int
    start::Int
    internal_deg::Int
end

@cache struct SERK2Cache{uType, rateType, uNoUnitsType, C <: SERK2ConstantCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Sᵢ::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    k::rateType
    constantcache::C
end

function alg_cache(
        alg::SERK2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    constantcache = SERK2ConstantCache(u)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    Sᵢ = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return SERK2Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::SERK2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SERK2ConstantCache(u)
end

mutable struct TSRKC3ConstantCache{zType, tTypeNoUnits} <: OrdinaryDiffEqConstantCache
    #to match the types to call maxeig!
    zprev::zType
    tsw0::tTypeNoUnits
    acoshtsw0::tTypeNoUnits
end
@cache struct TSRKC3Cache{uType, rateType, uNoUnitsType, C <: TSRKC3ConstantCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    gprev::uType
    gprev2::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    k::rateType
    constantcache::C
end

function alg_cache(
        alg::TSRKC3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tsw0 = tTypeNoUnits(1.25)
    acoshtsw0 = tTypeNoUnits(acosh(1.25))
    constantcache = TSRKC3ConstantCache(u, tsw0, acoshtsw0)
    gprev = zero(u)
    gprev2 = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return TSRKC3Cache(u, uprev, gprev, gprev2, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::TSRKC3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tsw0 = tTypeNoUnits(1.25)
    acoshtsw0 = tTypeNoUnits(acosh(1.25))
    return TSRKC3ConstantCache(u, tsw0, acoshtsw0)
end