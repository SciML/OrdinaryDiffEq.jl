struct PCEulerConstantCache <: StochasticDiffEqConstantCache end

function alg_cache(
        alg::PCEuler, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return PCEulerConstantCache()
end

@cache struct PCEulerCache{uType, rateType, rateNoiseType, rateNoiseCollectionType} <:
    StochasticDiffEqMutableCache
    utmp::uType
    ftmp::rateType
    gtmp::rateNoiseType
    gdWtmp::rateNoiseCollectionType
    bbprimetmp::rateType
end

function alg_cache(
        alg::PCEuler, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    utmp = zero(u)
    ftmp = zero(rate_prototype)
    gtmp = zero(noise_rate_prototype)
    bbprimetmp = zero(ftmp)
    if is_diagonal_noise(prob)
        gdWtmp = gtmp
    else
        gdWtmp = zero(rate_prototype)
    end
    return PCEulerCache(utmp, ftmp, gtmp, gdWtmp, bbprimetmp)
end
