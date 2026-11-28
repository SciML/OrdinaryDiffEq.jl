struct RandomAffineHeunConstantCache{randType, I1Type, I01Type, I2Type} <:
    StochasticDiffEqConstantCache
    wzero::randType
    wtmp::randType
    I1::I1Type
    I01::I01Type
    I2::I2Type
end
@cache struct RandomAffineHeunCache{uType, rateType, randType, I1Type, I01Type, I2Type} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    rtmp1::rateType
    rtmp2::rateType
    rtmp3::rateType
    wzero::randType
    wtmp::randType
    I1::I1Type
    I01::I01Type
    I2::I2Type
end

function affine_buffers(ΔW, dt)
    m = length(ΔW)
    w = zero(eltype(ΔW))
    return zeros(typeof(w * dt), m), zeros(typeof(w * dt^2), m),
        zeros(typeof(w^2 * dt^2), m, m)
end

function alg_cache(
        alg::RandomAffineHeun, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    I1, I01, I2 = affine_buffers(ΔW, dt)
    return RandomAffineHeunConstantCache(zero(ΔW), zero(ΔW), I1, I01, I2)
end

function alg_cache(
        alg::RandomAffineHeun, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    rtmp1 = zero(rate_prototype)
    rtmp2 = zero(rate_prototype)
    rtmp3 = zero(rate_prototype)
    I1, I01, I2 = affine_buffers(ΔW, dt)
    return RandomAffineHeunCache(
        u, uprev, tmp, rtmp1, rtmp2, rtmp3, zero(ΔW), zero(ΔW), I1, I01, I2
    )
end
