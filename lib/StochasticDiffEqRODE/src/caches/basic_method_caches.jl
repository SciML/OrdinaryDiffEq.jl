struct RandomEMConstantCache <: StochasticDiffEqConstantCache end
@cache struct RandomEMCache{uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    rtmp::rateType
end

function alg_cache(
        alg::RandomEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RandomEMConstantCache()
end

function alg_cache(
        alg::RandomEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    rtmp = zero(rate_prototype)
    return RandomEMCache(u, uprev, tmp, rtmp)
end

struct RandomTamedEMConstantCache <: StochasticDiffEqConstantCache end

@cache struct RandomTamedEMCache{uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    rtmp::rateType
end

function alg_cache(
        alg::RandomTamedEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RandomTamedEMConstantCache()
end

function alg_cache(
        alg::RandomTamedEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    rtmp = zero(rate_prototype)
    return RandomTamedEMCache(u, uprev, tmp, rtmp)
end

struct RandomHeunConstantCache <: StochasticDiffEqConstantCache end
@cache struct RandomHeunCache{uType, rateType, randType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    rtmp1::rateType
    rtmp2::rateType
    wtmp::randType
end

function alg_cache(
        alg::RandomHeun, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RandomHeunConstantCache()
end

function alg_cache(
        alg::RandomHeun, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    rtmp1 = zero(rate_prototype)
    rtmp2 = zero(rate_prototype)
    wtmp = zero(ΔW)
    return RandomHeunCache(u, uprev, tmp, rtmp1, rtmp2, wtmp)
end
