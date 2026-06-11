@cache mutable struct Tsit5Cache{uType, rateType, StageLimiter, StepLimiter, Thread, TmpC <: TmpCache} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    # Parameterized so the rate slots can be opted out (`precompute_initdt_cache`):
    # `TmpCache{uType, rateType, uNoUnitsType}` (rates held) or
    # `TmpCache{uType, Nothing, uNoUnitsType}` (rates skipped).
    tmp_cache::TmpC
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

@truncate_stacktrace Tsit5Cache 1

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
    # Tsit5 itself uses `tmp` (stage scratch), `tmp2` (embedded solution, was
    # `utilde`) and `atmp` (error-norm scaling). It needs no rate scratch of its
    # own — the rate buffers exist only so the initial-dt estimator can run
    # allocation-free. `precompute_initdt_cache` (default false) lets the user
    # opt into holding them; when false they are skipped and `initdt` allocates.
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits,
        alg.precompute_initdt_cache ? Val(true) : Val(false))
    return Tsit5Cache(
        u, uprev, k1, k2, k3, k4, k5, k6, k7, tmp_cache,
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
