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
