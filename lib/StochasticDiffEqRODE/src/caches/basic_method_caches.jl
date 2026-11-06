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

struct RandomTaylor15ConstantCache <: StochasticDiffEqConstantCache end
@cache struct RandomTaylor15Cache{uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    rtmp::rateType
    rtmp0::rateType
    rtmpp::rateType
    rtmpm::rateType
end

function warn_unresolved_grid(prob, t, dt)
    lo, hi = minmax(t, t + dt)
    if count(ti -> lo < ti < hi, prob.noise.t) == 0
        @warn "RandomTaylor15 is stepping a noise grid with no grid point inside the step, " *
            "so the step integrals reduce to the endpoint rule and the order drops to 1. " *
            "Use a finer noise grid or RandomEM."
    end
    return nothing
end

function alg_cache(
        alg::RandomTaylor15, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    warn_unresolved_grid(prob, t, dt)
    return RandomTaylor15ConstantCache()
end

function alg_cache(
        alg::RandomTaylor15, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    warn_unresolved_grid(prob, t, dt)
    tmp = zero(u)
    rtmp = zero(rate_prototype)
    rtmp0 = zero(rate_prototype)
    rtmpp = zero(rate_prototype)
    rtmpm = zero(rate_prototype)
    return RandomTaylor15Cache(u, uprev, tmp, rtmp, rtmp0, rtmpp, rtmpm)
end
