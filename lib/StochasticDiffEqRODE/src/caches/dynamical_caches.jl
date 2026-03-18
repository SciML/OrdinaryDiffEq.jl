mutable struct BAOABConstantCache{uType, uEltypeNoUnits} <: StochasticDiffEqConstantCache
    k::uType
    half::uEltypeNoUnits
    c1::uEltypeNoUnits
    c2::uEltypeNoUnits
end
@cache struct BAOABCache{uType, uEltypeNoUnits, rateNoiseType, uTypeCombined} <:
    StochasticDiffEqMutableCache
    utmp::uType
    dutmp::uType
    k::uType
    gtmp::uType
    noise::rateNoiseType
    half::uEltypeNoUnits
    c1::uEltypeNoUnits
    c2::uEltypeNoUnits
    tmp::uTypeCombined
end

function alg_cache(
        alg::BAOAB, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype.x[1])
    c1 = exp(-alg.gamma * dt)
    c2 = alg.scale_noise ? sqrt((1 - c1^2) / abs(dt)) : 1 # if scale_noise == false, c2 = 1
    return BAOABConstantCache(k, uEltypeNoUnits(1 // 2), uEltypeNoUnits(c1), uEltypeNoUnits(c2))
end

function alg_cache(
        alg::BAOAB, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dutmp = zero(u.x[1])
    utmp = zero(u.x[2])
    k = zero(rate_prototype.x[1])

    gtmp = zero(rate_prototype.x[1])
    noise = zero(rate_prototype.x[1])

    half = uEltypeNoUnits(1 // 2)
    c1 = exp(-alg.gamma * dt)
    c2 = alg.scale_noise ? sqrt((1 - c1^2) / abs(dt)) : 1 # if scale_noise == false, c2 = 1

    tmp = zero(u)

    return BAOABCache(
        utmp, dutmp, k, gtmp, noise, half, uEltypeNoUnits(c1), uEltypeNoUnits(c2), tmp
    )
end
