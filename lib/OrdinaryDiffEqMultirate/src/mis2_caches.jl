struct MIS2ConstantCache{TabType} <: OrdinaryDiffEqConstantCache
    tab::TabType
end

@cache mutable struct MIS2Cache{uType, rateType, uNoUnitsType, TabType} <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    atmp::uNoUnitsType
    v::uType
    offset::rateType
    k_fast::rateType
    Y::Vector{uType}
    fS::Vector{rateType}
    fsalfirst::rateType
    k::rateType
    tab::TabType
end

get_fsalfirstlast(cache::MIS2Cache, u) = (cache.fsalfirst, cache.k)

function alg_cache(
        alg::MIS2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = MIS2Tableau(eltype(u))
    s = length(tab.d)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    v = zero(u)
    offset = zero(rate_prototype)
    k_fast = zero(rate_prototype)
    Y = [zero(u) for _ in 1:s]
    fS = [zero(rate_prototype) for _ in 1:s]
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return MIS2Cache(u, uprev, tmp, atmp, v, offset, k_fast, Y, fS, fsalfirst, k, tab)
end

function alg_cache(
        alg::MIS2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MIS2ConstantCache(MIS2Tableau(eltype(u)))
end
