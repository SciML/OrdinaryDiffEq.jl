struct EMConstantCache <: StochasticDiffEqConstantCache end
@cache struct EMCache{uType, rateType, rateNoiseType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    rtmp1::rateType
    rtmp2::rateNoiseType
end

function alg_cache(
        alg::EM, prob, u, ΔW, ΔZ, p, rate_prototype, noise_rate_prototype,
        jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return EMConstantCache()
end

function alg_cache(
        alg::EM, prob, u, ΔW, ΔZ, p, rate_prototype, noise_rate_prototype,
        jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    rtmp1 = zero(rate_prototype)
    if noise_rate_prototype !== nothing
        rtmp2 = zero(noise_rate_prototype)
    else
        rtmp2 = nothing
    end
    return EMCache(u, uprev, tmp, rtmp1, rtmp2)
end

struct SplitEMConstantCache <: StochasticDiffEqConstantCache end
@cache struct SplitEMCache{uType, rateType, rateNoiseType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    rtmp1::rateType
    rtmp2::rateNoiseType
end

function alg_cache(
        alg::SplitEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SplitEMConstantCache()
end

function alg_cache(
        alg::SplitEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    rtmp1 = zero(rate_prototype)
    rtmp2 = zero(noise_rate_prototype)
    return SplitEMCache(u, uprev, tmp, rtmp1, rtmp2)
end

struct EulerHeunConstantCache <: StochasticDiffEqConstantCache end
@cache struct EulerHeunCache{uType, rateType, rateNoiseType, rateNoiseCollectionType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    ftmp1::rateType
    ftmp2::rateType
    nrtmp::rateNoiseCollectionType
    gtmp1::rateNoiseType
    gtmp2::rateNoiseType
end

function alg_cache(
        alg::EulerHeun, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return EulerHeunConstantCache()
end

function alg_cache(
        alg::EulerHeun, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    ftmp1 = zero(rate_prototype)
    ftmp2 = zero(rate_prototype)
    nrtmp = zero(rate_prototype)
    gtmp1 = zero(noise_rate_prototype)
    gtmp2 = zero(noise_rate_prototype)
    return EulerHeunCache(u, uprev, tmp, ftmp1, ftmp2, nrtmp, gtmp1, gtmp2)
end

struct SimplifiedEMConstantCache <: StochasticDiffEqConstantCache end
@cache struct SimplifiedEMCache{randType, uType, rateType, rateNoiseType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    _dW::randType
    rtmp1::rateType
    rtmp2::rateNoiseType
end

function alg_cache(
        alg::SimplifiedEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SimplifiedEMConstantCache()
end

function alg_cache(
        alg::SimplifiedEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
    else
        _dW = zero(ΔW)
    end

    rtmp1 = zero(rate_prototype)
    rtmp2 = zero(noise_rate_prototype)
    return SimplifiedEMCache(u, uprev, _dW, rtmp1, rtmp2)
end

struct RKMilConstantCache <: StochasticDiffEqConstantCache end
@cache struct RKMilCache{uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    du1::rateType
    du2::rateType
    K::rateType
    tmp::uType
    L::rateType
end

function alg_cache(
        alg::RKMil, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RKMilConstantCache()
end

function alg_cache(
        alg::RKMil, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    K = zero(rate_prototype)
    tmp = zero(u)
    L = zero(rate_prototype)
    return RKMilCache(u, uprev, du1, du2, K, tmp, L)
end

struct RKMilCommuteConstantCache{JalgType} <: StochasticDiffEqConstantCache
    Jalg::JalgType
end
@cache struct RKMilCommuteCache{uType, rateType, rateNoiseType, JalgType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    du1::rateType
    du2::rateType
    K::rateType
    gtmp::rateNoiseType
    L::rateNoiseType
    Jalg::JalgType
    mil_correction::rateType
    Kj::uType
    Dgj::rateNoiseType
    tmp::uType
end

function alg_cache(
        alg::RKMilCommute, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    Jalg = get_Jalg(ΔW, dt, prob, alg)
    return RKMilCommuteConstantCache{typeof(Jalg)}(Jalg)
end

function alg_cache(
        alg::RKMilCommute, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    K = zero(rate_prototype)
    gtmp = zero(noise_rate_prototype)
    L = zero(noise_rate_prototype)
    tmp = zero(rate_prototype)
    Jalg = get_Jalg(ΔW, dt, prob, alg)
    mil_correction = zero(rate_prototype)
    Kj = zero(u)
    Dgj = zero(noise_rate_prototype)
    return RKMilCommuteCache(u, uprev, du1, du2, K, gtmp, L, Jalg, mil_correction, Kj, Dgj, tmp)
end
