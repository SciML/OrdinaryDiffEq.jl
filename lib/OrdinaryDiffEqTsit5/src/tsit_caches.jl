@cache struct Tsit5Cache{
        uType, rateType, uNoUnitsType, StageLimiter, StepLimiter,
        Thread,
    } <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

"""Concrete cache for Tsit5 with `Vector{Float64}` and default limiters (0 type parameters)."""
struct Tsit5CacheVF64 <: OrdinaryDiffEqMutableCache
    u::Vector{Float64}
    uprev::Vector{Float64}
    k1::Vector{Float64}
    k2::Vector{Float64}
    k3::Vector{Float64}
    k4::Vector{Float64}
    k5::Vector{Float64}
    k6::Vector{Float64}
    k7::Vector{Float64}
    utilde::Vector{Float64}
    tmp::Vector{Float64}
    atmp::Vector{Float64}
    stage_limiter!::typeof(trivial_limiter!)
    step_limiter!::typeof(trivial_limiter!)
    thread::False
end

full_cache(c::Tsit5CacheVF64) = (c.u, c.uprev, c.k1, c.k2, c.k3, c.k4, c.k5, c.k6, c.k7, c.utilde, c.tmp, c.atmp)

const Tsit5CacheType = Union{Tsit5Cache, Tsit5CacheVF64}

function alg_cache(
        alg::Tsit5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    if u isa Vector{Float64} && rate_prototype isa Vector{Float64} &&
            uEltypeNoUnits === Float64 &&
            alg.stage_limiter! isa typeof(trivial_limiter!) &&
            alg.step_limiter! isa typeof(trivial_limiter!) &&
            alg.thread isa False
        return Tsit5CacheVF64(
            u, uprev, k1, k2, k3, k4, k5, k6, k7, utilde, tmp, atmp,
            alg.stage_limiter!, alg.step_limiter!, alg.thread
        )
    end
    return Tsit5Cache(
        u, uprev, k1, k2, k3, k4, k5, k6, k7, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread
    )
end

get_fsalfirstlast(cache::Tsit5CacheType, u) = (cache.k1, cache.k7)

function alg_cache(
        alg::Tsit5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return Tsit5ConstantCache()
end
