struct MRABConstantCache{TabType} <: OrdinaryDiffEqConstantCache
    tab::TabType
end

@cache mutable struct MRABCache{uType, rateType, uNoUnitsType, TabType} <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    atmp::uNoUnitsType
    k_slow::rateType
    k_fast::rateType
    F_history::Vector{rateType}
    fsalfirst::rateType
    k::rateType
    tab::TabType
end

get_fsalfirstlast(cache::MRABCache, u) = (cache.fsalfirst, cache.k)

function alg_cache(
        alg::MRAB, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    k_slow = zero(rate_prototype)
    k_fast = zero(rate_prototype)
    F_history = [zero(rate_prototype) for _ in 1:(alg.k)]
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    tab = MRABTableau(alg.k, eltype(u))
    return MRABCache(u, uprev, tmp, atmp, k_slow, k_fast, F_history, fsalfirst, k, tab)
end

function alg_cache(
        alg::MRAB, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MRABConstantCache(MRABTableau(alg.k, eltype(u)))
end
