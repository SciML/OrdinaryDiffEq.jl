@mutable_cache mutable struct Tsit5Cache{
        uType, rateType, uNoUnitsType, StageLimiter, StepLimiter,
        Thread,
    } <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
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

# Buffer swap is beneficial for arrays larger than ~1000 elements where memcpy
# of 7 k-arrays per step dominates. For smaller arrays, allocation overhead
# of 7 fresh arrays per step outweighs the copy savings.
const _K_SWAP_THRESHOLD = 1000
OrdinaryDiffEqCore.supports_k_swap(cache::Tsit5Cache) = length(cache.k1) >= _K_SWAP_THRESHOLD

function OrdinaryDiffEqCore.swap_k_buffers!(integrator, cache::Tsit5Cache)
    # Capture FSAL reference before allocating new arrays: old k7 has the
    # correct FSAL data from the just-completed step. By pointing fsallast
    # to it, update_fsal! will naturally copy it into the new k1 (fsalfirst)
    # via its existing recursivecopy!(fsalfirst, fsallast). No extra copy needed.
    integrator.fsallast = cache.k7

    # Allocate fresh k arrays for the next step
    cache.k1 = similar(cache.k7)
    cache.k2 = similar(cache.k7)
    cache.k3 = similar(cache.k7)
    cache.k4 = similar(cache.k7)
    cache.k5 = similar(cache.k7)
    cache.k6 = similar(cache.k7)
    cache.k7 = similar(cache.k7)

    # Update integrator.k references to new arrays
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    integrator.k[6] = cache.k6
    integrator.k[7] = cache.k7

    # fsalfirst points to new k1 (update_fsal! will copy fsallast into it)
    integrator.fsalfirst = cache.k1

    return nothing
end

function alg_cache(
        alg::Tsit5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return Tsit5ConstantCache()
end
