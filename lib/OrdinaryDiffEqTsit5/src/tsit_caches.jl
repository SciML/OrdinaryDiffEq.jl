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
    # Parameterized so the rate slots can be opted in (`preallocate_initdt_buffers`):
    # `TmpCache{uType, Nothing, uNoUnitsType, Nothing}` (default, rates skipped) or
    # `TmpCache{uType, rateType, uNoUnitsType, Nothing}` (rates held).
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
    # `utilde`) and `atmp` (error-norm scaling) — migrated fields, so the array
    # count matches the historical cache exactly. It has no safe rate donors
    # (`k1..k7` all feed the interpolant), so the initdt rate buffers are
    # allocated only when the user opts in via `preallocate_initdt_buffers`
    # (default false); otherwise those slots are `nothing` and `initdt`
    # allocates its rate temporaries at call time.
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits,
        alg.preallocate_initdt_buffers ? Val(true) : Val(false))
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
