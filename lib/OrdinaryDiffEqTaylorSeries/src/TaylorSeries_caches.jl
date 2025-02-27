@cache struct ExplicitTaylor2Cache{
    uType, rateType, uNoUnitsType, StageLimiter, StepLimiter,
    Thread} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::ExplicitTaylor2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    ExplicitTaylor2Cache(u, uprev, k1, k2, k3, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread)
end
struct ExplicitTaylor2ConstantCache <: OrdinaryDiffEqConstantCache end
function alg_cache(alg::ExplicitTaylor2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    ExplicitTaylor2ConstantCache()
end
# FSAL currently not used, providing dummy implementation to satisfy the interface
get_fsalfirstlast(cache::ExplicitTaylor2Cache, u) = (cache.k1, cache.k1)

@cache struct ExplicitTaylorCache{
    P, uType, rateType, StageLimiter, StepLimiter,
    Thread} <: OrdinaryDiffEqMutableCache
    order::Val{P}
    u::uType
    uprev::uType
    us::NTuple{P, uType}
    ks::NTuple{P, rateType}
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::ExplicitTaylor{P}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {P, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # ks: normalized derivatives for f, starting from 0
    ks = ntuple(i -> zero(rate_prototype), Val(P))
    # us: normalized derivatives for u, starting from 1
    us = ntuple(i -> zero(u), Val(P))
    ExplicitTaylorCache(Val(P), u, uprev, us, ks, alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

struct ExplicitTaylorConstantCache{P} <: OrdinaryDiffEqConstantCache end
function alg_cache(alg::ExplicitTaylor{P}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {P, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    ExplicitTaylorConstantCache{P}()
end

# FSAL currently not used, providing dummy implementation to satisfy the interface
get_fsalfirstlast(cache::ExplicitTaylorCache, u) = (cache.ks[1], cache.ks[1])


# Differential Algebriac Equation Taylor Series
@cache mutable struct DAETSCache{
    uType, rateType, uNoUnitsType, tTypeNoUnits, StageLimiter, StepLimiter,
    Thread, MatType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    utilde::uType
    tmp::uType
    atmp::uType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
    xcur::uType
    xTS::uType
    xtrial::uType
    htrial::tTypeNoUnits
    e::tTypeNoUnits
    Σ::MatType        # Signature matrix
    c::Vector{Int}    # Offsets for equations
    d::Vector{Int}    # Offsets for variables
    J::MatType        # Jacobian
end

function alg_cache(alg::DAETS, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck, ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    xcur = similar(u)
    xTS = similar(u)
    xtrial = similar(u)
    htrial = dt
    e = zero(tTypeNoUnits)
    tmp = similar(u)
    atmp = similar(u)
    utilde = similar(u)

    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)  

    fill!(xcur, 0.0)
    fill!(xTS, 0.0)
    fill!(xtrial, 0.0)
    fill!(tmp, 0.0)
    fill!(atmp, 0.0)
    fill!(utilde, 0.0)

    # DAETS specific fields
    n = length(u)
    Σ = zeros(n, n)  # Empty signature matrix
    c = zeros(Int, n)  # Equation offsets
    d = zeros(Int, n)  # Variable offsets
    J = zeros(n, n)  # Jacobian

    # Return cache
    return DAETSCache{
        typeof(u), typeof(rate_prototype), typeof(atmp), tTypeNoUnits,
        typeof(alg.stage_limiter!), typeof(alg.step_limiter!), typeof(alg.thread),
        typeof(Σ)
    }(
        u, uprev, k1, k2, k3, utilde, tmp, atmp, 
        alg.stage_limiter!, alg.step_limiter!, alg.thread,
        xcur, xTS, xtrial, htrial, e,
        Σ, c, d, J 
    )
end
get_fsalfirstlast(cache::DAETSCache, u) = (cache.k1, cache.k1)