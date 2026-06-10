abstract type StabilizedRKMutableCache <: OrdinaryDiffEqMutableCache end
get_fsalfirstlast(cache::StabilizedRKMutableCache, u) = (cache.fsalfirst, cache.k)

mutable struct ROCK2ConstantCache{T, T2, zType} <: OrdinaryDiffEqConstantCache
    ms::NTuple{46, Int}
    fp1::NTuple{46, T}
    fp2::NTuple{46, T}
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
    ms::NTuple{50, Int}
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

mutable struct RKMC2ConstantCache{zType, T} <: OrdinaryDiffEqConstantCache
    zprev::zType
    mdeg::Int
    min_stage::Int
    max_stage::Int
    w0::T
    w1::T
end

@cache struct RKMC2Cache{uType, rateType, uNoUnitsType, C <: RKMC2ConstantCache} <:
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
        alg::RKMC2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    w = zero(tTypeNoUnits)
    constantcache = RKMC2ConstantCache(zero(u), 0, alg.min_stages, alg.max_stages, w, w)
    gprev = zero(u)
    gprev2 = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return RKMC2Cache(u, uprev, gprev, gprev2, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::RKMC2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    w = zero(tTypeNoUnits)
    return RKMC2ConstantCache(zero(u), 0, alg.min_stages, alg.max_stages, w, w)
end

mutable struct ESERK4ConstantCache{T, zType} <: OrdinaryDiffEqConstantCache
    ms::NTuple{46, Int}
    Cᵤ::NTuple{4, Int}
    Cₑ::NTuple{4, Int}
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
    constantcache = ESERK4ConstantCache(constvalue(uBottomEltypeNoUnits), u)
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
    return ESERK4ConstantCache(constvalue(uBottomEltypeNoUnits), u)
end

mutable struct ESERK5ConstantCache{T, zType} <: OrdinaryDiffEqConstantCache
    ms::NTuple{49, Int}
    Cᵤ::NTuple{5, Int}
    Cₑ::NTuple{5, Int}
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
    constantcache = ESERK5ConstantCache(constvalue(uBottomEltypeNoUnits), u)
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
    return ESERK5ConstantCache(constvalue(uBottomEltypeNoUnits), u)
end

mutable struct SERK2ConstantCache{T, zType} <: OrdinaryDiffEqConstantCache
    ms::NTuple{11, Int}
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
    constantcache = SERK2ConstantCache(constvalue(uBottomEltypeNoUnits), u)
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
    return SERK2ConstantCache(constvalue(uBottomEltypeNoUnits), u)
end

mutable struct TSRKC2ConstantCache{zType} <: OrdinaryDiffEqConstantCache
    #to match the types to call maxeig!
    zprev::zType
end
@cache struct TSRKC2Cache{uType, rateType, uNoUnitsType, C <: TSRKC2ConstantCache} <:
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
        alg::TSRKC2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    constantcache = TSRKC2ConstantCache(u)
    gprev = zero(u)
    gprev2 = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return TSRKC2Cache(u, uprev, gprev, gprev2, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::TSRKC2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return TSRKC2ConstantCache(u)
end

mutable struct TSRKC3ConstantCache{zType} <: OrdinaryDiffEqConstantCache
    #to match the types to call maxeig!
    zprev::zType
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
    constantcache = TSRKC3ConstantCache(u)
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
    return TSRKC3ConstantCache(u)
end

mutable struct RKL1ConstantCache{zType} <: OrdinaryDiffEqConstantCache
    zprev::zType
    mdeg::Int
    min_stage::Int
    max_stage::Int
end

function RKL1ConstantCache(u, min_stage, max_stage)
    return RKL1ConstantCache(zero(u), 3, min_stage, max_stage)
end

@cache struct RKL1Cache{uType, rateType, uNoUnitsType, C <: RKL1ConstantCache} <:
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
        alg::RKL1, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck, ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    min_stage, max_stage = _rkl_clamp_odd_stages(alg.min_stages, alg.max_stages)
    constantcache = RKL1ConstantCache(u, min_stage, max_stage)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return RKL1Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::RKL1, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck, ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    min_stage, max_stage = _rkl_clamp_odd_stages(alg.min_stages, alg.max_stages)
    return RKL1ConstantCache(u, min_stage, max_stage)
end

mutable struct RKL2ConstantCache{zType} <: OrdinaryDiffEqConstantCache
    zprev::zType
    mdeg::Int
    min_stage::Int
    max_stage::Int
end

function RKL2ConstantCache(u, min_stage, max_stage)
    return RKL2ConstantCache(zero(u), 3, min_stage, max_stage)
end

@cache struct RKL2Cache{uType, rateType, uNoUnitsType, C <: RKL2ConstantCache} <:
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
        alg::RKL2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck, ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    min_stage, max_stage = _rkl_clamp_odd_stages(alg.min_stages, alg.max_stages)
    constantcache = RKL2ConstantCache(u, min_stage, max_stage)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return RKL2Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::RKL2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck, ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    min_stage, max_stage = _rkl_clamp_odd_stages(alg.min_stages, alg.max_stages)
    return RKL2ConstantCache(u, min_stage, max_stage)
end

mutable struct RKG1ConstantCache{zType} <: OrdinaryDiffEqConstantCache
    zprev::zType
    mdeg::Int
    min_stage::Int
    max_stage::Int
end

function RKG1ConstantCache(u, min_stage, max_stage)
    return RKG1ConstantCache(zero(u), 2, min_stage, max_stage)
end

@cache struct RKG1Cache{uType, rateType, uNoUnitsType, C <: RKG1ConstantCache} <:
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
        alg::RKG1, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck, ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    min_stage = max(2, alg.min_stages)
    max_stage = max(alg.max_stages, min_stage)
    constantcache = RKG1ConstantCache(u, min_stage, max_stage)
    uᵢ₋₁ = zero(u); uᵢ₋₂ = zero(u); tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits); recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype); k = zero(rate_prototype)
    return RKG1Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::RKG1, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck, ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    min_stage = max(2, alg.min_stages)
    max_stage = max(alg.max_stages, min_stage)
    return RKG1ConstantCache(u, min_stage, max_stage)
end

mutable struct RKG2ConstantCache{zType} <: OrdinaryDiffEqConstantCache
    zprev::zType
    mdeg::Int
    min_stage::Int
    max_stage::Int
end

function RKG2ConstantCache(u, min_stage, max_stage)
    return RKG2ConstantCache(zero(u), 3, min_stage, max_stage)
end

@cache struct RKG2Cache{uType, rateType, uNoUnitsType, C <: RKG2ConstantCache} <:
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
        alg::RKG2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck, ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    min_stage = max(3, alg.min_stages)
    max_stage = max(alg.max_stages, min_stage)
    constantcache = RKG2ConstantCache(u, min_stage, max_stage)
    uᵢ₋₁ = zero(u); uᵢ₋₂ = zero(u); tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits); recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype); k = zero(rate_prototype)
    return RKG2Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::RKG2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck, ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    min_stage = max(3, alg.min_stages)
    max_stage = max(alg.max_stages, min_stage)
    return RKG2ConstantCache(u, min_stage, max_stage)
end
