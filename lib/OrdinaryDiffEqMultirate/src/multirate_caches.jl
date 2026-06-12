struct MREEFConstantCache{T} <: OrdinaryDiffEqConstantCache
    T::T  # pre-allocated extrapolation table: Vector of length `order`
end

@cache mutable struct MREEFCache{uType, rateType, uNoUnitsType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    atmp::uNoUnitsType
    k_slow::rateType
    k_fast::rateType
    T::Array{uType, 1}
    fsalfirst::rateType
    k::rateType
end

get_fsalfirstlast(cache::MREEFCache, u) = (cache.fsalfirst, cache.k)

function alg_cache(
        alg::MREEF, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    k_slow = zero(rate_prototype)
    k_fast = zero(rate_prototype)
    T = [zero(u) for _ in 1:(alg.order)]
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return MREEFCache(u, uprev, tmp, atmp, k_slow, k_fast, T, fsalfirst, k)
end

function alg_cache(
        alg::MREEF, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    T = Vector{typeof(u)}(undef, alg.order)
    return MREEFConstantCache(T)
end

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

struct MISConstantCache{TabType} <: OrdinaryDiffEqConstantCache
    tab::TabType
end

@cache mutable struct MISCache{uType, rateType, uNoUnitsType, TabType} <:
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

get_fsalfirstlast(cache::MISCache, u) = (cache.fsalfirst, cache.k)

function alg_cache(
        alg::MIS, u, rate_prototype, ::Type{uEltypeNoUnits},
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
    return MISCache(u, uprev, tmp, atmp, v, offset, k_fast, Y, fS, fsalfirst, k, tab)
end

function alg_cache(
        alg::MIS, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MISConstantCache(MIS2Tableau(eltype(u)))
end
