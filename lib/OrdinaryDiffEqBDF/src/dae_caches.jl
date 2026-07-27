abstract type DAEBDFMutableCache <: OrdinaryDiffEqMutableCache end
function get_fsalfirstlast(cache::DAEBDFMutableCache, u)
    return (cache.fsalfirst, du_alias_or_new(cache.nlsolver, cache.fsalfirst))
end

@cache mutable struct DImplicitEulerCache{uType, rateType, N, TmpC <: TmpCache} <:
    DAEBDFMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    # Unified scratch: only the cache-level `atmp` migrated here (the Newton
    # buffers live on the nlsolver and are off limits; `k₁`/`k₂` feed
    # `integrator.k` and are never legal donors). `tmp`/`tmp2`/`weight` are
    # `nothing`; the rate slots stay opted out since DImplicitEuler's
    # positional constructor exposes no `preallocate_initdt_buffers` knob:
    # `TmpCache{Nothing, Nothing, uNoUnitsType, Nothing}`.
    tmp_cache::TmpC
    k₁::rateType
    k₂::rateType
    nlsolver::N
end

# Not FSAL
get_fsalfirstlast(cache::DImplicitEulerCache, u) = (nothing, nothing)

mutable struct DImplicitEulerConstantCache{N} <: OrdinaryDiffEqConstantCache
    nlsolver::N
end

function alg_cache(
        alg::DImplicitEuler, du, u, res_prototype, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1, 1
    α = 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, res_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, α, Val(false), verbose
    )

    return DImplicitEulerConstantCache(nlsolver)
end

function alg_cache(
        alg::DImplicitEuler, du, u, res_prototype, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1, 1
    α = 1
    k₁ = zero(rate_prototype)
    k₂ = zero(rate_prototype)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, res_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, α, Val(true), verbose
    )

    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    # `atmp` is the only migrated scratch field; everything else stays
    # `nothing`, so the array count matches the historical cache exactly.
    tmp_cache = TmpCache(nothing, nothing, atmp, nothing, nothing, nothing)

    return DImplicitEulerCache(u, uprev, uprev2, tmp_cache, k₁, k₂, nlsolver)
end

@cache mutable struct DABDF2ConstantCache{N, dtType, rate_prototype} <:
    OrdinaryDiffEqConstantCache
    nlsolver::N
    eulercache::DImplicitEulerConstantCache
    dtₙ₋₁::dtType
    fsalfirstprev::rate_prototype
end

function alg_cache(
        alg::DABDF2, du, u, res_prototype, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = Int64(1) // 1, 1
    α = Int64(1) // 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, res_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, α, Val(false), verbose
    )
    eulercache = DImplicitEulerConstantCache(nlsolver)

    dtₙ₋₁ = one(dt)
    fsalfirstprev = rate_prototype

    return DABDF2ConstantCache(nlsolver, eulercache, dtₙ₋₁, fsalfirstprev)
end

@cache mutable struct DABDF2Cache{uType, rateType, N, dtType, TmpC <: TmpCache} <:
    DAEBDFMutableCache
    uₙ::uType
    uₙ₋₁::uType
    uₙ₋₂::uType
    fsalfirst::rateType
    fsalfirstprev::rateType
    # Unified scratch: only the cache-level `atmp` migrated here (shared with
    # the bootstrap eulercache, preserving the historical aliasing). The other
    # slots stay `nothing`; the rate slots stay opted out since DABDF2's
    # positional constructor exposes no `preallocate_initdt_buffers` knob:
    # `TmpCache{Nothing, Nothing, uNoUnitsType, Nothing}`.
    tmp_cache::TmpC
    nlsolver::N
    eulercache::DImplicitEulerCache
    dtₙ₋₁::dtType
end

function alg_cache(
        alg::DABDF2, du, u, res_prototype, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = Int64(1) // 1, 1
    α = Int64(1) // 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, res_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, α, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    fsalfirstprev = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    k₁ = zero(rate_prototype)
    k₂ = zero(rate_prototype)

    # One TmpCache shared by the outer cache and the bootstrap eulercache —
    # both wrap the same `atmp` array (preserving the historical aliasing), so
    # the array count matches the historical cache exactly.
    tmp_cache = TmpCache(nothing, nothing, atmp, nothing, nothing, nothing)

    eulercache = DImplicitEulerCache(u, uprev, uprev2, tmp_cache, k₁, k₂, nlsolver)

    dtₙ₋₁ = one(dt)

    return DABDF2Cache(
        u, uprev, uprev2, fsalfirst, fsalfirstprev, tmp_cache,
        nlsolver, eulercache, dtₙ₋₁
    )
end

@cache mutable struct DFBDFConstantCache{
        MO, N, tsType, tType, uType, uuType, coeffType,
        EEstType, rType, wType, fdWeightsType,
    } <:
    OrdinaryDiffEqConstantCache
    nlsolver::N
    ts::tsType
    ts_tmp::tsType
    t_old::tType
    u_history::uuType
    order::Int
    prev_order::Int
    u₀::uType
    u_corrector::uuType
    bdf_coeffs::coeffType
    max_order::Val{MO}
    nconsteps::Int
    consfailcnt::Int
    qwait::Int # countdown to next order change consideration (CVODE-style)
    terkm2::EEstType
    terkm1::EEstType
    terk::EEstType
    terkp1::EEstType
    r::rType
    weights::wType
    iters_from_event::Int
    fd_weights::fdWeightsType
end

function alg_cache(
        alg::DFBDF{MO}, du, u, res_prototype, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits,
        uprev, uprev2, f, t, dt, reltol, p, calck, ::Val{false}, verbose
    ) where {MO}
    γ, c = 1.0, 1.0
    max_order = MO
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    bdf_coeffs = _make_bdf_coeffs_dfbdf()
    ts = zero(Vector{typeof(t)}(undef, max_order + 2)) #ts is the successful past points, it will be updated after successful step
    ts_tmp = similar(ts)

    if u isa Number
        u_history = zeros(eltype(u), max_order + 2)
        u_corrector = zeros(eltype(u), max_order + 2)
    else
        u_history = [zero(u) for _ in 1:(max_order + 2)]
        u_corrector = [zero(u) for _ in 1:(max_order + 2)]
    end
    order = 1
    prev_order = 1
    terkm2 = tTypeNoUnits(1)
    terkm1 = tTypeNoUnits(1)
    terk = tTypeNoUnits(1)
    terkp1 = tTypeNoUnits(1)
    r = zero(Vector{typeof(t)}(undef, max_order + 2))
    weights = zero(Vector{typeof(t)}(undef, max_order + 2))
    weights[1] = 1
    nconsteps = 0
    consfailcnt = 0
    qwait = 3 # order + 2, matching nconsteps >= order + 2 for failure-free runs
    t_old = zero(t)
    iters_from_event = 0
    u₀ = zero(u)

    fd_weights = zeros(typeof(t), max_order + 1, max_order + 1)

    return DFBDFConstantCache(
        nlsolver, ts, ts_tmp, t_old, u_history, order, prev_order, u₀,
        u_corrector, bdf_coeffs, Val(MO), nconsteps, consfailcnt, qwait, terkm2,
        terkm1, terk, terkp1, r, weights, iters_from_event, fd_weights
    )
end

@cache mutable struct DFBDFCache{
        MO, N, rateType, tsType, tType, uType,
        uuType, coeffType, EEstType, rType, wType, fdWeightsType,
        TmpC <: TmpCache,
    } <:
    DAEBDFMutableCache
    fsalfirst::rateType
    nlsolver::N
    ts::tsType
    ts_tmp::tsType
    t_old::tType
    u_history::uuType
    order::Int
    prev_order::Int
    u_corrector::uuType
    u₀::uType
    bdf_coeffs::coeffType
    max_order::Val{MO}
    nconsteps::Int
    consfailcnt::Int
    qwait::Int # countdown to next order change consideration (CVODE-style)
    # Unified scratch: the former inline `tmp` (corrector RHS, write-first
    # every step) lives in the `tmp` slot and `atmp` in `atmp`; no `utilde`
    # existed so `tmp2` stays `nothing`. The error-estimate buffers
    # `terk_tmp`/`terkp1_tmp` keep their dedicated fields. Rate slots stay
    # opted out (no `preallocate_initdt_buffers` knob on DFBDF):
    # `TmpCache{uType, Nothing, uNoUnitsType, Nothing}`.
    tmp_cache::TmpC
    terkm2::EEstType
    terkm1::EEstType
    terk::EEstType
    terkp1::EEstType
    terk_tmp::uType
    terkp1_tmp::uType
    r::rType
    weights::wType
    equi_ts::tsType
    iters_from_event::Int
    dense::Vector{uType}
    fd_weights::fdWeightsType
end

function alg_cache(
        alg::DFBDF{MO}, du, u, res_prototype, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {MO}
    γ, c = 1.0, 1.0
    fsalfirst = zero(rate_prototype)
    max_order = MO
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    bdf_coeffs = _make_bdf_coeffs_dfbdf()
    ts = Vector{typeof(t)}(undef, max_order + 2) #ts is the successful past points, it will be updated after successful step
    u_history = [zero(u) for _ in 1:(max_order + 2)]
    order = 1
    prev_order = 1
    u_corrector = [zero(u) for _ in 1:(max_order + 2)]
    recursivefill!(ts, zero(t))
    terkm2 = tTypeNoUnits(1)
    terkm1 = tTypeNoUnits(1)
    terk = tTypeNoUnits(1)
    terkp1 = tTypeNoUnits(1)
    terk_tmp = similar(u)
    terkp1_tmp = similar(u)
    r = Vector{typeof(t)}(undef, max_order + 2)
    weights = Vector{typeof(t)}(undef, max_order + 2)
    recursivefill!(r, zero(t))
    recursivefill!(weights, zero(t))
    weights[1] = 1
    nconsteps = 0
    consfailcnt = 0
    qwait = 3 # order + 2, matching nconsteps >= order + 2 for failure-free runs
    t_old = zero(t)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, zero(uEltypeNoUnits))
    u₀ = similar(u)
    equi_ts = similar(ts)
    tmp = similar(u)
    ts_tmp = similar(ts)
    iters_from_event = 0

    dense = [zero(u) for _ in 1:(2 * (max_order + 1))]  # first half for integrator.k, second half as scratch

    fd_weights = zeros(typeof(t), max_order + 1, max_order + 1)

    # Migrated fields only (`tmp` → tmp, `atmp` → atmp) — the array count
    # matches the historical cache exactly.
    tmp_cache = TmpCache(tmp, nothing, atmp, nothing, nothing, nothing)

    return DFBDFCache(
        fsalfirst, nlsolver, ts, ts_tmp, t_old, u_history, order, prev_order,
        u_corrector, u₀, bdf_coeffs, Val(MO), nconsteps, consfailcnt, qwait, tmp_cache,
        terkm2, terkm1, terk, terkp1, terk_tmp, terkp1_tmp, r, weights, equi_ts,
        iters_from_event, dense, fd_weights
    )
end
