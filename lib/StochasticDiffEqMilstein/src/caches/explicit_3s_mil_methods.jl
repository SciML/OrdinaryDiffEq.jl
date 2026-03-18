struct WangLi3SMil_AConstantCache <: StochasticDiffEqConstantCache end
@cache struct WangLi3SMil_ACache{uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    k::rateType
    k1::rateType
    tmp::uType
end

function alg_cache(
        alg::WangLi3SMil_A, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return WangLi3SMil_AConstantCache()
end

function alg_cache(
        alg::WangLi3SMil_A, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    k1 = zero(rate_prototype)
    tmp = zero(u)
    return WangLi3SMil_ACache(u, uprev, k, k1, tmp)
end

struct WangLi3SMil_BConstantCache <: StochasticDiffEqConstantCache end
@cache struct WangLi3SMil_BCache{uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    k::rateType
    k1::rateType
    tmp::uType
end

function alg_cache(
        alg::WangLi3SMil_B, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return WangLi3SMil_BConstantCache()
end

function alg_cache(
        alg::WangLi3SMil_B, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    k1 = zero(rate_prototype)
    tmp = zero(u)
    return WangLi3SMil_BCache(u, uprev, k, k1, tmp)
end

struct WangLi3SMil_CConstantCache <: StochasticDiffEqConstantCache end
@cache struct WangLi3SMil_CCache{uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    k::rateType
    k1::rateType
    tmp::uType
end

function alg_cache(
        alg::WangLi3SMil_C, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return WangLi3SMil_CConstantCache()
end

function alg_cache(
        alg::WangLi3SMil_C, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    k1 = zero(rate_prototype)
    tmp = zero(u)
    return WangLi3SMil_CCache(u, uprev, k, k1, tmp)
end

struct WangLi3SMil_DConstantCache <: StochasticDiffEqConstantCache end
@cache struct WangLi3SMil_DCache{uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    k::rateType
    k1::rateType
    tmp::uType
end

function alg_cache(
        alg::WangLi3SMil_D, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return WangLi3SMil_DConstantCache()
end

function alg_cache(
        alg::WangLi3SMil_D, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    k1 = zero(rate_prototype)
    tmp = zero(u)
    return WangLi3SMil_DCache(u, uprev, k, k1, tmp)
end

struct WangLi3SMil_EConstantCache <: StochasticDiffEqConstantCache end
@cache struct WangLi3SMil_ECache{uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    k::rateType
    k1::rateType
    tmp::uType
end

function alg_cache(
        alg::WangLi3SMil_E, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return WangLi3SMil_EConstantCache()
end

function alg_cache(
        alg::WangLi3SMil_E, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    k1 = zero(rate_prototype)
    tmp = zero(u)
    return WangLi3SMil_ECache(u, uprev, k, k1, tmp)
end

struct WangLi3SMil_FConstantCache <: StochasticDiffEqConstantCache end
@cache struct WangLi3SMil_FCache{uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    k::rateType
    k1::rateType
    tmp::uType
end

function alg_cache(
        alg::WangLi3SMil_F, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return WangLi3SMil_FConstantCache()
end

function alg_cache(
        alg::WangLi3SMil_F, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    k1 = zero(rate_prototype)
    tmp = zero(u)
    return WangLi3SMil_FCache(u, uprev, k, k1, tmp)
end
