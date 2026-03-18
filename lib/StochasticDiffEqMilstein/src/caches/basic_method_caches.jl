struct RKMilGeneralConstantCache{JalgType} <: StochasticDiffEqConstantCache
    Jalg::JalgType
end

@cache struct RKMilGeneralCache{uType, rateType, rateNoiseType, JalgType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    du1::rateType
    du2::rateType
    K::uType
    L::rateNoiseType
    mil_correction::uType
    ggprime::rateNoiseType
    Jalg::JalgType
end

function alg_cache(
        alg::RKMilGeneral, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    Jalg = get_Jalg(DeltaW, dt, prob, alg)
    return RKMilGeneralConstantCache{typeof(Jalg)}(Jalg)
end

function alg_cache(
        alg::RKMilGeneral, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    K = zero(u)
    L = zero(noise_rate_prototype)
    mil_correction = zero(u)
    ggprime = zero(noise_rate_prototype)
    Jalg = get_Jalg(DeltaW, dt, prob, alg)
    return RKMilGeneralCache{
        typeof(u), typeof(rate_prototype), typeof(noise_rate_prototype), typeof(Jalg),
    }(
        u, uprev, tmp, du1, du2, K, L, mil_correction, ggprime, Jalg
    )
end
