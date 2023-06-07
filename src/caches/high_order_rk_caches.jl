@cache struct TanYam7Cache{uType, rateType, uNoUnitsType, TabType, StageLimiter,
    StepLimiter, Thread} <:
              OrdinaryDiffEqMutableCache
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
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    k::rateType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::TanYam7, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = TanYam7ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
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
    utilde = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    k = zero(rate_prototype)
    TanYam7Cache(u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, utilde, tmp, atmp, k,
        tab, alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

function alg_cache(alg::TanYam7, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    TanYam7ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct DP8Cache{uType, rateType, uNoUnitsType, TabType, StageLimiter, StepLimiter,
    Thread} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
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
    kupdate::rateType
    udiff::rateType
    bspl::rateType
    dense_tmp3::rateType
    dense_tmp4::rateType
    dense_tmp5::rateType
    dense_tmp6::rateType
    dense_tmp7::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::DP8, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
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
    kupdate = zero(rate_prototype)
    utilde = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    k13 = zero(rate_prototype)
    k14 = zero(rate_prototype)
    k15 = zero(rate_prototype)
    k16 = zero(rate_prototype)
    udiff = zero(rate_prototype)
    bspl = zero(rate_prototype)
    # dense_tmp1 = udiff
    # dense_tmp2 = bspl
    dense_tmp3 = zero(rate_prototype)
    dense_tmp4 = zero(rate_prototype)
    dense_tmp5 = zero(rate_prototype)
    dense_tmp6 = zero(rate_prototype)
    dense_tmp7 = zero(rate_prototype)
    tab = DP8ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    DP8Cache(u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15,
        k16, kupdate,
        udiff, bspl, dense_tmp3, dense_tmp4, dense_tmp5, dense_tmp6, dense_tmp7,
        utilde, tmp, atmp, tab, alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

function alg_cache(alg::DP8, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    DP8ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct TsitPap8Cache{uType, rateType, uNoUnitsType, TabType, StageLimiter,
    StepLimiter, Thread} <:
              OrdinaryDiffEqMutableCache
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
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    k::rateType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::TsitPap8, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = TsitPap8ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
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
    utilde = zero(u)
    k = zero(rate_prototype)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    TsitPap8Cache(u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, utilde,
        tmp, atmp, k, tab, alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

function alg_cache(alg::TsitPap8, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    TsitPap8ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct PFRK87Cache{uType, rateType, uNoUnitsType, TabType, StageLimiter, StepLimiter,
    Thread} <:
              OrdinaryDiffEqMutableCache
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
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    k::rateType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::PFRK87, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = PFRK87ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
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
    utilde = zero(u)
    k = zero(rate_prototype)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    PFRK87Cache(u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, utilde,
        tmp, atmp, k, tab, alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

function alg_cache(alg::PFRK87, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    PFRK87ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end
