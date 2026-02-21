@cache struct Vern6Cache{
        uType, rateType, uNoUnitsType, TabType, StageLimiter, StepLimiter,
        Thread, L,
    } <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    k9::rateType
    utilde::uType
    tmp::uType
    rtmp::rateType
    atmp::uNoUnitsType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
    lazy::L
end

get_fsalfirstlast(cache::Vern6Cache, u) = (cache.k1, cache.k9)

function alg_cache(
        alg::Vern6, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Vern6Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = k2
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = k3
    k9 = zero(rate_prototype)
    utilde = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    rtmp = uEltypeNoUnits === eltype(u) ? utilde : zero(rate_prototype)
    return Vern6Cache(
        u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, k9, utilde, tmp, rtmp, atmp, tab,
        alg.stage_limiter!, alg.step_limiter!, alg.thread, alg.lazy
    )
end

struct Vern6ConstantCache{TabType, L} <: OrdinaryDiffEqConstantCache
    tab::TabType
    lazy::L
end

function alg_cache(
        alg::Vern6, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Vern6Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return Vern6ConstantCache(tab, alg.lazy)
end

@cache struct Vern7Cache{
        uType, rateType, uNoUnitsType, StageLimiter, StepLimiter,
        Thread, L,
    } <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    k9::rateType
    k10::rateType
    utilde::uType
    tmp::uType
    rtmp::rateType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
    lazy::L
end

"""Concrete cache for Vern7 with `Vector{Float64}` and default limiters (0 type parameters)."""
struct Vern7CacheVF64 <: OrdinaryDiffEqMutableCache
    u::Vector{Float64}
    uprev::Vector{Float64}
    k1::Vector{Float64}
    k2::Vector{Float64}
    k3::Vector{Float64}
    k4::Vector{Float64}
    k5::Vector{Float64}
    k6::Vector{Float64}
    k7::Vector{Float64}
    k8::Vector{Float64}
    k9::Vector{Float64}
    k10::Vector{Float64}
    utilde::Vector{Float64}
    tmp::Vector{Float64}
    rtmp::Vector{Float64}
    atmp::Vector{Float64}
    stage_limiter!::typeof(trivial_limiter!)
    step_limiter!::typeof(trivial_limiter!)
    thread::False
    lazy::Val{true}
end

full_cache(c::Vern7CacheVF64) = (c.u, c.uprev, c.k1, c.k2, c.k3, c.k4, c.k5, c.k6, c.k7, c.k8, c.k9, c.k10, c.utilde, c.tmp, c.rtmp, c.atmp)

const Vern7CacheType = Union{Vern7Cache, Vern7CacheVF64}

# fake values since non-FSAL method
get_fsalfirstlast(cache::Vern7CacheType, u) = (nothing, nothing)

function alg_cache(
        alg::Vern7, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = k2
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    k9 = zero(rate_prototype)
    k10 = k2
    utilde = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    rtmp = uEltypeNoUnits === eltype(u) ? utilde : zero(rate_prototype)
    if u isa Vector{Float64} && rate_prototype isa Vector{Float64} &&
            uEltypeNoUnits === Float64 &&
            alg.stage_limiter! isa typeof(trivial_limiter!) &&
            alg.step_limiter! isa typeof(trivial_limiter!) &&
            alg.thread isa False && alg.lazy isa Val{true}
        return Vern7CacheVF64(
            u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, utilde, tmp, rtmp, atmp,
            alg.stage_limiter!, alg.step_limiter!, alg.thread, alg.lazy
        )
    end
    return Vern7Cache(
        u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, utilde, tmp, rtmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread, alg.lazy
    )
end

struct Vern7ConstantCache{L} <: OrdinaryDiffEqConstantCache
    lazy::L
end

function alg_cache(
        alg::Vern7, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return Vern7ConstantCache(alg.lazy)
end

@cache struct Vern8Cache{
        uType, rateType, uNoUnitsType, TabType, StageLimiter, StepLimiter,
        Thread, L,
    } <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    k9::rateType
    k10::rateType
    k11::rateType
    k12::rateType
    k13::rateType
    utilde::uType
    tmp::uType
    rtmp::rateType
    atmp::uNoUnitsType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
    lazy::L
end

# fake values since non-FSAL method
get_fsalfirstlast(cache::Vern8Cache, u) = (nothing, nothing)

function alg_cache(
        alg::Vern8, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Vern8Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = k2
    k4 = zero(rate_prototype)
    k5 = k2
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    tmp = zero(u)
    k9 = zero(rate_prototype)
    k10 = zero(rate_prototype)
    k11 = zero(rate_prototype)
    k12 = zero(rate_prototype)
    k13 = k4
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    rtmp = uEltypeNoUnits === eltype(u) ? utilde : zero(rate_prototype)
    return Vern8Cache(
        u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, utilde,
        tmp, rtmp, atmp, tab, alg.stage_limiter!, alg.step_limiter!, alg.thread, alg.lazy
    )
end

struct Vern8ConstantCache{TabType, L} <: OrdinaryDiffEqConstantCache
    tab::TabType
    lazy::L
end

function alg_cache(
        alg::Vern8, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Vern8Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return Vern8ConstantCache(tab, alg.lazy)
end

@cache struct Vern9Cache{
        uType, rateType, uNoUnitsType, StageLimiter, StepLimiter,
        Thread, L,
    } <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    k9::rateType
    k10::rateType
    k11::rateType
    k12::rateType
    k13::rateType
    k14::rateType
    k15::rateType
    k16::rateType
    utilde::uType
    tmp::uType
    rtmp::rateType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
    lazy::L
end

# fake values since non-FSAL method
get_fsalfirstlast(cache::Vern9Cache, u) = (nothing, nothing)

function alg_cache(
        alg::Vern9, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = k2
    k4 = zero(rate_prototype)
    k5 = k3
    k6 = zero(rate_prototype)
    k7 = k4
    k8 = k5
    k9 = zero(rate_prototype)
    k10 = zero(rate_prototype)
    k11 = zero(rate_prototype)
    k12 = zero(rate_prototype)
    k13 = zero(rate_prototype)
    k14 = zero(rate_prototype)
    k15 = zero(rate_prototype)
    k16 = k6
    utilde = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    rtmp = uEltypeNoUnits === eltype(u) ? utilde : zero(rate_prototype)
    return Vern9Cache(
        u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15,
        k16, utilde, tmp, rtmp, atmp, alg.stage_limiter!, alg.step_limiter!,
        alg.thread, alg.lazy
    )
end

struct Vern9ConstantCache{L} <: OrdinaryDiffEqConstantCache
    lazy::L
end

function alg_cache(
        alg::Vern9, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return Vern9ConstantCache(alg.lazy)
end

@cache struct RKV76IIaCache{
        uType, rateType, uNoUnitsType, TabType, StageLimiter, StepLimiter,
        Thread, L,
    } <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    k9::rateType
    k10::rateType
    utilde::uType
    tmp::uType
    rtmp::rateType
    atmp::uNoUnitsType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
    lazy::L
end

# fake values since non-FSAL method
get_fsalfirstlast(cache::RKV76IIaCache, u) = (nothing, nothing)

function alg_cache(
        alg::RKV76IIa, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = RKV76IIaTableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = k2
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = k3
    k9 = zero(rate_prototype)
    k10 = k4
    utilde = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    rtmp = uEltypeNoUnits === eltype(u) ? utilde : zero(rate_prototype)
    return RKV76IIaCache(
        u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, utilde, tmp, rtmp, atmp, tab,
        alg.stage_limiter!, alg.step_limiter!, alg.thread, alg.lazy
    )
end

struct RKV76IIaConstantCache{TabType, L} <: OrdinaryDiffEqConstantCache
    tab::TabType
    lazy::L
end

function alg_cache(
        alg::RKV76IIa, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = RKV76IIaTableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return RKV76IIaConstantCache(tab, alg.lazy)
end
