struct LambaEMConstantCache <: StochasticDiffEqConstantCache end
@cache struct LambaEMCache{uType, rateType, rateNoiseType, randType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    du1::rateType
    du2::rateType
    K::rateType
    tmp::uType
    L::rateNoiseType
    gtmp::rateNoiseType
    dW_cache::randType
end

function alg_cache(
        alg::LambaEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return LambaEMConstantCache()
end

function alg_cache(
        alg::LambaEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    K = zero(rate_prototype)
    tmp = zero(u)
    L = zero(noise_rate_prototype)
    gtmp = zero(noise_rate_prototype)
    if is_diagonal_noise(prob)
        dW_cache = nothing
    else
        dW_cache = zero(ΔW)
    end
    return LambaEMCache(u, uprev, du1, du2, K, tmp, L, gtmp, dW_cache)
end

struct LambaEulerHeunConstantCache <: StochasticDiffEqConstantCache end
@cache struct LambaEulerHeunCache{uType, rateType, rateNoiseType, randType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    du1::rateType
    du2::rateType
    K::rateType
    tmp::uType
    L::rateNoiseType
    gtmp::rateNoiseType
    dW_cache::randType
end

function alg_cache(
        alg::LambaEulerHeun, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return LambaEulerHeunConstantCache()
end

function alg_cache(
        alg::LambaEulerHeun, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    K = zero(rate_prototype)
    tmp = zero(u)
    L = zero(noise_rate_prototype)
    gtmp = zero(noise_rate_prototype)
    if is_diagonal_noise(prob)
        dW_cache = nothing
    else
        dW_cache = zero(ΔW)
    end
    return LambaEulerHeunCache(u, uprev, du1, du2, K, tmp, L, gtmp, dW_cache)
end
