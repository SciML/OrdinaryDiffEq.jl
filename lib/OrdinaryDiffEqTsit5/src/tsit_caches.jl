using ConcreteStructs

@concrete struct Tsit5Cache <: OrdinaryDiffEqMutableCache
    u
    uprev
    k1
    k2
    k3
    k4
    k5
    k6
    k7
    utilde
    tmp
    atmp
    stage_limiter!
    step_limiter!
    thread
end

function alg_cache(
        alg::Tsit5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return Tsit5Cache(
        u, uprev, k1, k2, k3, k4, k5, k6, k7, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread
    )
end

get_fsalfirstlast(cache::Tsit5Cache, u) = (cache.k1, cache.k7)

function alg_cache(
        alg::Tsit5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return Tsit5ConstantCache()
end
