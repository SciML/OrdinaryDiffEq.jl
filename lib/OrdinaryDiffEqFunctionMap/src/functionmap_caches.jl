@cache struct FunctionMapCache{uType, rateType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::rateType
end

function alg_cache(alg::FunctionMap, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    FunctionMapCache(u, uprev,
        FunctionMap_scale_by_time(alg) ? rate_prototype :
        (eltype(u) <: Enum ? copy(u) : zero(u)))
end

struct FunctionMapConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::FunctionMap, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    FunctionMapConstantCache()
end

isdiscretecache(cache::Union{FunctionMapCache, FunctionMapConstantCache}) = true
