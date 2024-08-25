@cache struct Tsit5Cache{uType, rateType, uNoUnitsType, StageLimiter, StepLimiter,
    Thread} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k::Vector{rateType}
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::Tsit5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = [zero(rate_prototype) for _ in 1:7]
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    Tsit5Cache(u, uprev, k, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

get_fsalfirstlast(cache::Tsit5Cache, u) = (cache.k[1], cache.k[7])

function alg_cache(alg::Tsit5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    Tsit5ConstantCache()
end
