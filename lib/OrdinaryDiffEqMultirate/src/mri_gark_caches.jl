struct MRIGARKERK22ConstantCache{TabType} <: OrdinaryDiffEqConstantCache
    tab::TabType
end

@cache mutable struct MRIGARKERK22Cache{uType, rateType, uNoUnitsType, TabType} <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    atmp::uNoUnitsType
    v::uType
    f1eval::rateType
    fS_yn::rateType
    fS_Y2::rateType
    Y2::uType
    fsalfirst::rateType
    k::rateType
    tab::TabType
end

get_fsalfirstlast(cache::MRIGARKERK22Cache, u) = (cache.fsalfirst, cache.k)

function alg_cache(
        alg::MRIGARKERK22a, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    v = zero(u)
    f1eval = zero(rate_prototype)
    fS_yn = zero(rate_prototype)
    fS_Y2 = zero(rate_prototype)
    Y2 = zero(u)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    tab = MRIGARKERK22aTableau(eltype(u))
    return MRIGARKERK22Cache(
        u, uprev, tmp, atmp, v, f1eval, fS_yn, fS_Y2, Y2, fsalfirst, k, tab
    )
end

function alg_cache(
        alg::MRIGARKERK22a, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MRIGARKERK22ConstantCache(MRIGARKERK22aTableau(eltype(u)))
end

function alg_cache(
        alg::MRIGARKERK22b, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    v = zero(u)
    f1eval = zero(rate_prototype)
    fS_yn = zero(rate_prototype)
    fS_Y2 = zero(rate_prototype)
    Y2 = zero(u)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    tab = MRIGARKERK22bTableau(eltype(u))
    return MRIGARKERK22Cache(
        u, uprev, tmp, atmp, v, f1eval, fS_yn, fS_Y2, Y2, fsalfirst, k, tab
    )
end

function alg_cache(
        alg::MRIGARKERK22b, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MRIGARKERK22ConstantCache(MRIGARKERK22bTableau(eltype(u)))
end
