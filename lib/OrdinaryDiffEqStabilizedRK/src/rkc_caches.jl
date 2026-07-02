abstract type StabilizedRKMutableCache <: OrdinaryDiffEqMutableCache end
get_fsalfirstlast(cache::StabilizedRKMutableCache, u) = (cache.fsalfirst, cache.k)

# All StabilizedRK mutable caches share the same scratch layout: an inline `tmp`
# (state scratch, also reused as the eigenvector buffer by `maxeig!`) and `atmp`
# (error-norm scratch), with no `utilde` — so the `tmp2` and `weight` slots are
# opted out (`nothing`) and the raw `TmpCache` constructor is used (array count
# identical to the historical caches). The fsal buffers `fsalfirst`/`k` are
# exposed through `integrator.k` (dense output), so there are no legal rate
# donors; the initdt rate slots are freshly allocated only when the algorithm
# opts in via the `preallocate_initdt_buffers` trait (none of the StabilizedRK
# algorithms carry that field today, so the trait folds to `false` and the
# slots are `nothing`, letting `initdt` allocate at call time instead).
function stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return _stabilized_rk_tmp_cache(
        tmp, atmp, rate_prototype, Val(preallocate_initdt_buffers(alg))
    )
end
function _stabilized_rk_tmp_cache(
        tmp, atmp, rate_prototype, ::Val{need_rates}
    ) where {need_rates}
    # `need_rates` is a type parameter, so these ternaries fold at compile time.
    rate_tmp = need_rates ? zero(rate_prototype) : nothing
    rate_tmp2 = need_rates ? zero(rate_prototype) : nothing
    return TmpCache(tmp, nothing, atmp, nothing, rate_tmp, rate_tmp2)
end

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
@cache struct ROCK2Cache{uType, rateType, C <: ROCK2ConstantCache, TmpC <: TmpCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    tmp_cache::TmpC
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
    tmp_cache = stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return ROCK2Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, tmp_cache, fsalfirst, k, constantcache)
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

@cache struct ROCK4Cache{uType, rateType, C <: ROCK4ConstantCache, TmpC <: TmpCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    uᵢ₋₃::uType
    tmp_cache::TmpC
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
    tmp_cache = stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return ROCK4Cache(
        u, uprev, uᵢ₋₁, uᵢ₋₂, uᵢ₋₃, tmp_cache, fsalfirst, k, constantcache
    )
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
@cache struct RKCCache{uType, rateType, C <: RKCConstantCache, TmpC <: TmpCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    gprev::uType
    gprev2::uType
    tmp_cache::TmpC
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
    tmp_cache = stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return RKCCache(u, uprev, gprev, gprev2, tmp_cache, fsalfirst, k, constantcache)
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

@cache struct RKMC2Cache{uType, rateType, C <: RKMC2ConstantCache, TmpC <: TmpCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    gprev::uType
    gprev2::uType
    tmp_cache::TmpC
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
    tmp_cache = stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return RKMC2Cache(u, uprev, gprev, gprev2, tmp_cache, fsalfirst, k, constantcache)
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

@cache struct ESERK4Cache{uType, rateType, C <: ESERK4ConstantCache, TmpC <: TmpCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Sᵢ::uType
    tmp_cache::TmpC
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
    tmp_cache = stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return ESERK4Cache(
        u, uprev, uᵢ, uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp_cache, fsalfirst, k, constantcache
    )
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

@cache struct ESERK5Cache{uType, rateType, C <: ESERK5ConstantCache, TmpC <: TmpCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Sᵢ::uType
    tmp_cache::TmpC
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
    tmp_cache = stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return ESERK5Cache(
        u, uprev, uᵢ, uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp_cache, fsalfirst, k, constantcache
    )
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

@cache struct SERK2Cache{uType, rateType, C <: SERK2ConstantCache, TmpC <: TmpCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Sᵢ::uType
    tmp_cache::TmpC
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
    tmp_cache = stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return SERK2Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp_cache, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::SERK2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SERK2ConstantCache(constvalue(uBottomEltypeNoUnits), u)
end

mutable struct TSRKC2ConstantCache{zType, tTypeNoUnits} <: OrdinaryDiffEqConstantCache
    #to match the types to call maxeig!
    zprev::zType
    tsw0::tTypeNoUnits
    acoshtsw0::tTypeNoUnits
end
@cache struct TSRKC2Cache{uType, rateType, C <: TSRKC2ConstantCache, TmpC <: TmpCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    gprev::uType
    gprev2::uType
    tmp_cache::TmpC
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
    tsw0 = tTypeNoUnits(1.1)
    acoshtsw0 = tTypeNoUnits(acosh(1.1))
    constantcache = TSRKC2ConstantCache(u, tsw0, acoshtsw0)
    gprev = zero(u)
    gprev2 = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    tmp_cache = stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return TSRKC2Cache(u, uprev, gprev, gprev2, tmp_cache, fsalfirst, k, constantcache)
end

function alg_cache(
        alg::TSRKC2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tsw0 = tTypeNoUnits(1.1)
    acoshtsw0 = tTypeNoUnits(acosh(1.1))
    return TSRKC2ConstantCache(u, tsw0, acoshtsw0)
end

mutable struct TSRKC3ConstantCache{zType, tTypeNoUnits} <: OrdinaryDiffEqConstantCache
    #to match the types to call maxeig!
    zprev::zType
    tsw0::tTypeNoUnits
    acoshtsw0::tTypeNoUnits
end
@cache struct TSRKC3Cache{uType, rateType, C <: TSRKC3ConstantCache, TmpC <: TmpCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    gprev::uType
    gprev2::uType
    tmp_cache::TmpC
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
    tmp_cache = stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return TSRKC3Cache(u, uprev, gprev, gprev2, tmp_cache, fsalfirst, k, constantcache)
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

mutable struct RKL1ConstantCache{zType} <: OrdinaryDiffEqConstantCache
    zprev::zType
    mdeg::Int
    min_stage::Int
    max_stage::Int
end

function RKL1ConstantCache(u, min_stage, max_stage)
    return RKL1ConstantCache(zero(u), 3, min_stage, max_stage)
end

@cache struct RKL1Cache{uType, rateType, C <: RKL1ConstantCache, TmpC <: TmpCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    tmp_cache::TmpC
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
    tmp_cache = stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return RKL1Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, tmp_cache, fsalfirst, k, constantcache)
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

@cache struct RKL2Cache{uType, rateType, C <: RKL2ConstantCache, TmpC <: TmpCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    tmp_cache::TmpC
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
    tmp_cache = stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return RKL2Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, tmp_cache, fsalfirst, k, constantcache)
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

@cache struct RKG1Cache{uType, rateType, C <: RKG1ConstantCache, TmpC <: TmpCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    tmp_cache::TmpC
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
    tmp_cache = stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return RKG1Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, tmp_cache, fsalfirst, k, constantcache)
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

@cache struct RKG2Cache{uType, rateType, C <: RKG2ConstantCache, TmpC <: TmpCache} <:
    StabilizedRKMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    tmp_cache::TmpC
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
    tmp_cache = stabilized_rk_tmp_cache(alg, tmp, atmp, rate_prototype)
    return RKG2Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, tmp_cache, fsalfirst, k, constantcache)
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
