struct MRIGARKConstantCache{TabType} <: OrdinaryDiffEqConstantCache
    tab::TabType
end

@cache mutable struct MRIGARKCache{uType, rateType, uNoUnitsType, TabType} <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    atmp::uNoUnitsType
    v::uType
    f1eval::rateType
    slow_f::rateType
    Y::Vector{uType}
    fS::Vector{rateType}
    fsalfirst::rateType
    k::rateType
    tab::TabType
end

get_fsalfirstlast(cache::MRIGARKCache, u) = (cache.fsalfirst, cache.k)

function alg_cache(
        alg::Union{MRIGARKERK22a, MRIGARKERK22b}, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = mri_gark_tableau(alg, eltype(u))
    s = length(tab.Γ)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    v = zero(u)
    f1eval = zero(rate_prototype)
    slow_f = zero(rate_prototype)
    Y = [zero(u) for _ in 1:s]
    fS = [zero(rate_prototype) for _ in 1:s]
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return MRIGARKCache(u, uprev, tmp, atmp, v, f1eval, slow_f, Y, fS, fsalfirst, k, tab)
end

function alg_cache(
        alg::Union{MRIGARKERK22a, MRIGARKERK22b}, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MRIGARKConstantCache(mri_gark_tableau(alg, eltype(u)))
end
