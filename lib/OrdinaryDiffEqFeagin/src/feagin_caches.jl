abstract type FeaginCache <: OrdinaryDiffEqMutableCache end
get_fsalfirstlast(cache::FeaginCache, u) = (cache.fsalfirst, cache.k)

@cache struct Feagin10Cache{uType, uNoUnitsType, rateType, TabType, StepLimiter} <:
    FeaginCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    k9::rateType
    k10::rateType
    k11::rateType
    k12::rateType
    k13::rateType
    k14::rateType
    k15::rateType
    k16::rateType
    k17::rateType
    tmp::uType
    atmp::uNoUnitsType
    k::rateType
    tab::TabType
    step_limiter!::StepLimiter
end

function alg_cache(
        alg::Feagin10, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Feagin10ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    k9 = zero(rate_prototype)
    k10 = zero(rate_prototype)
    k11 = zero(rate_prototype)
    k12 = zero(rate_prototype)
    k13 = zero(rate_prototype)
    k14 = zero(rate_prototype)
    k15 = zero(rate_prototype)
    k16 = zero(rate_prototype)
    k17 = zero(rate_prototype)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    k = zero(rate_prototype)

    return Feagin10Cache(
        u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14,
        k15, k16, k17, tmp, atmp, k, tab, alg.step_limiter!
    )
end

function alg_cache(
        alg::Feagin10, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return Feagin10ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct Feagin12Cache{uType, uNoUnitsType, rateType, TabType, StepLimiter} <:
    FeaginCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    k9::rateType
    k10::rateType
    k11::rateType
    k12::rateType
    k13::rateType
    k14::rateType
    k15::rateType
    k16::rateType
    k17::rateType
    k18::rateType
    k19::rateType
    k20::rateType
    k21::rateType
    k22::rateType
    k23::rateType
    k24::rateType
    k25::rateType
    tmp::uType
    atmp::uNoUnitsType
    k::rateType
    tab::TabType
    step_limiter!::StepLimiter
end

function alg_cache(
        alg::Feagin12, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Feagin12ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    k9 = zero(rate_prototype)
    k10 = zero(rate_prototype)
    k11 = zero(rate_prototype)
    k12 = zero(rate_prototype)
    k13 = zero(rate_prototype)
    k14 = zero(rate_prototype)
    k15 = zero(rate_prototype)
    k16 = zero(rate_prototype)
    k17 = zero(rate_prototype)
    k18 = zero(rate_prototype)
    k19 = zero(rate_prototype)
    k20 = zero(rate_prototype)
    k21 = zero(rate_prototype)
    k22 = zero(rate_prototype)
    k23 = zero(rate_prototype)
    k24 = zero(rate_prototype)
    k25 = zero(rate_prototype)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    k = zero(rate_prototype)

    return Feagin12Cache(
        u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14,
        k15, k16, k17, k18, k19, k20, k21, k22, k23, k24,
        k25, tmp, atmp, k, tab, alg.step_limiter!
    )
end

function alg_cache(
        alg::Feagin12, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return Feagin12ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct Feagin14Cache{uType, uNoUnitsType, rateType, TabType, StepLimiter} <:
    FeaginCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    k9::rateType
    k10::rateType
    k11::rateType
    k12::rateType
    k13::rateType
    k14::rateType
    k15::rateType
    k16::rateType
    k17::rateType
    k18::rateType
    k19::rateType
    k20::rateType
    k21::rateType
    k22::rateType
    k23::rateType
    k24::rateType
    k25::rateType
    k26::rateType
    k27::rateType
    k28::rateType
    k29::rateType
    k30::rateType
    k31::rateType
    k32::rateType
    k33::rateType
    k34::rateType
    k35::rateType
    tmp::uType
    atmp::uNoUnitsType
    k::rateType
    tab::TabType
    step_limiter!::StepLimiter
end

function alg_cache(
        alg::Feagin14, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Feagin14ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    k9 = zero(rate_prototype)
    k10 = zero(rate_prototype)
    k11 = zero(rate_prototype)
    k12 = zero(rate_prototype)
    k13 = zero(rate_prototype)
    k14 = zero(rate_prototype)
    k15 = zero(rate_prototype)
    k16 = zero(rate_prototype)
    k17 = zero(rate_prototype)
    k18 = zero(rate_prototype)
    k19 = zero(rate_prototype)
    k20 = zero(rate_prototype)
    k21 = zero(rate_prototype)
    k22 = zero(rate_prototype)
    k23 = zero(rate_prototype)
    k24 = zero(rate_prototype)
    k25 = zero(rate_prototype)
    k26 = zero(rate_prototype)
    k27 = zero(rate_prototype)
    k28 = zero(rate_prototype)
    k29 = zero(rate_prototype)
    k30 = zero(rate_prototype)
    k31 = zero(rate_prototype)
    k32 = zero(rate_prototype)
    k33 = zero(rate_prototype)
    k34 = zero(rate_prototype)
    k35 = zero(rate_prototype)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    k = zero(rate_prototype)

    return Feagin14Cache(
        u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14,
        k15, k16,
        k17, k18, k19, k20, k21, k22, k23, k24, k25, k26, k27, k28, k29, k30,
        k31, k32, k33, k34, k35, tmp, atmp, k, tab, alg.step_limiter!
    )
end

function alg_cache(
        alg::Feagin14, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return Feagin14ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end
