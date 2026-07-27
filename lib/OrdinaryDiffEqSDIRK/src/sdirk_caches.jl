abstract type SDIRKMutableCache <: OrdinaryDiffEqMutableCache end
abstract type SDIRKConstantCache <: OrdinaryDiffEqConstantCache end
function get_fsalfirstlast(cache::SDIRKMutableCache, u)
    return (cache.fsalfirst, du_alias_or_new(cache.nlsolver, cache.fsalfirst))
end

# Pure SDIRK types (non-IMEX) unified into ESDIRKIMEXCache
const _PureSDIRKAlg = Union{
    OrdinaryDiffEqNewtonNonAdaptiveSDIRKAlgorithm,
    OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm,
    ImplicitEuler, Trapezoid, CFNLIRK3,
}

# step_limiter! accessor — only some pure SDIRK algorithms have the field
_esdirk_step_limiter!(alg::OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm) = alg.step_limiter!
_esdirk_step_limiter!(alg::Union{ImplicitMidpoint, SDIRK2, TRBDF2, ImplicitEuler, Trapezoid}) = alg.step_limiter!
_esdirk_step_limiter!(alg) = trivial_limiter!

# smooth_est accessor — only adaptive algorithms carry this flag
_esdirk_smooth_est(alg::OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm) = alg.smooth_est
_esdirk_smooth_est(alg::OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm) = alg.smooth_est
_esdirk_smooth_est(alg) = false

"""
    ESDIRKIMEXConstantCache <: OrdinaryDiffEqConstantCache

Out-of-place solver cache for the ESDIRK-IMEX methods. Holds the nonlinear solver
`nlsolver`, the Butcher tableau `tab`, and the extra history slots `uprev3` /
`tprev2` used by the embedded error estimate. Declared public so downstream IMEX
solvers can reuse the ESDIRK-IMEX step.
"""
mutable struct ESDIRKIMEXConstantCache{Tab, N, U3, T2} <: OrdinaryDiffEqConstantCache
    nlsolver::N
    tab::Tab
    uprev3::U3
    tprev2::T2
end

function ESDIRKIMEXConstantCache(nlsolver, tab)
    return ESDIRKIMEXConstantCache(nlsolver, tab, nothing, nothing)
end

"""
    ESDIRKIMEXCache <: SDIRKMutableCache

In-place solver cache for the ESDIRK-IMEX methods. Holds the stage values `zs`,
stage-derivative buffers `ks`, error temporary `atmp`, nonlinear solver
`nlsolver`, tableau `tab`, step limiter, and the extra history slots
(`uprev2`/`uprev3`/`tprev2`) and `algebraic_vars` mask.
"""
mutable struct ESDIRKIMEXCache{
        uType, rateType, N, Tab, kType, StepLimiter, U2, AV, U3, T2, TmpC <: TmpCache,
    } <: SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    zs::Vector{uType}
    ks::Vector{kType}
    # Unified scratch: only the cache-level `atmp` migrated here (Newton caches
    # carry no cache-level state scratch — the Newton buffers live on the
    # nlsolver and are off limits), so `tmp`/`tmp2`/`weight` are `nothing`.
    # The rate slots can be opted in via `preallocate_initdt_buffers`:
    # `TmpCache{Nothing, Nothing, uNoUnitsType, Nothing}` (default) or
    # `TmpCache{Nothing, rateType, uNoUnitsType, Nothing}` (rates held).
    tmp_cache::TmpC
    nlsolver::N
    tab::Tab
    step_limiter!::StepLimiter
    uprev2::U2
    algebraic_vars::AV
    uprev3::U3
    tprev2::T2
end

# Compatibility shim: external constructors (OrdinaryDiffEqBDF's ABDF2
# first-step bootstrap builds this cache positionally) still pass a bare
# `atmp` array in the slot that now holds the unified scratch struct. Wrap it
# so the same array keeps serving as the error-norm scratch.
function ESDIRKIMEXCache(
        u, uprev, fsalfirst, zs, ks, atmp::AbstractArray, nlsolver, tab,
        step_limiter!, uprev2, algebraic_vars, uprev3, tprev2
    )
    tmp_cache = TmpCache(nothing, nothing, atmp, nothing, nothing, nothing)
    return ESDIRKIMEXCache(
        u, uprev, fsalfirst, zs, ks, tmp_cache, nlsolver, tab, step_limiter!,
        uprev2, algebraic_vars, uprev3, tprev2
    )
end

function ESDIRKIMEXCache(u, uprev, fsalfirst, zs, ks, tmp_cache, nlsolver, tab, step_limiter!)
    return ESDIRKIMEXCache(
        u, uprev, fsalfirst, zs, ks, tmp_cache, nlsolver, tab, step_limiter!,
        nothing, nothing, nothing, nothing
    )
end

function full_cache(c::ESDIRKIMEXCache)
    # Opted-out TmpCache slots appear as `nothing`; full_cache consumers
    # (resize! & friends) skip those.
    base = (
        c.u, c.uprev, c.fsalfirst, c.zs...,
        OrdinaryDiffEqCore.tmp_cache_buffers(c.tmp_cache)...,
    )
    if eltype(c.ks) !== Nothing
        base = tuple(base..., c.ks...)
    end
    if c.uprev2 !== nothing
        base = tuple(base..., c.uprev2)
    end
    if c.uprev3 !== nothing
        base = tuple(base..., c.uprev3)
    end
    return base
end

function OrdinaryDiffEqCore.strip_cache(cache::ESDIRKIMEXCache)
    s = length(cache.zs)
    return ConstructionBase.constructorof(typeof(cache))(
        nothing, nothing, nothing,
        Vector{Nothing}(undef, s),
        Vector{Nothing}(undef, s),
        TmpCache(nothing, nothing, nothing, nothing, nothing, nothing),
        nothing, nothing, nothing, nothing, nothing, nothing, nothing
    )
end

# Rate-scratch buffers for the initial-dt estimator, gated at the type level on
# the algorithm's `preallocate_initdt_buffers` trait (always `false` for the
# current SDIRK algorithms, which have no such field — the slots stay
# `nothing` and `initdt` allocates its rate temporaries at call time). There is
# no safe rate donor to alias: `fsalfirst` feeds the Hermite interpolant.
_initdt_rate_buffers(rate_prototype, ::Val{true}) =
    (zero(rate_prototype), zero(rate_prototype))
_initdt_rate_buffers(rate_prototype, ::Val{false}) = (nothing, nothing)

# `atmp` is the only cache-level scratch this cache ever had, so it is the only
# migrated slot — `tmp`/`tmp2` stay `nothing` rather than allocating state
# scratch the historical cache didn't have (the Newton state on the nlsolver is
# deliberately not aliased here). Net array count is unchanged from master.
function _esdirk_tmp_cache(alg, u, rate_prototype, ::Type{uEltypeNoUnits}) where {uEltypeNoUnits}
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    rate_tmp, rate_tmp2 = _initdt_rate_buffers(
        rate_prototype, Val(OrdinaryDiffEqCore.preallocate_initdt_buffers(alg))
    )
    return TmpCache(nothing, nothing, atmp, nothing, rate_tmp, rate_tmp2)
end

function alg_cache(
        alg::OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = ESDIRKIMEXTableau(alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ = tab.Ai[tab.s, tab.s]
    c = tab.nlsolver_init_c
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    return ESDIRKIMEXConstantCache(nlsolver, tab)
end

function alg_cache(
        alg::OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = ESDIRKIMEXTableau(alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ = tab.Ai[tab.s, tab.s]
    c = tab.nlsolver_init_c
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    s = tab.s
    if f isa SplitFunction
        ks = [zero(u) for _ in 1:s]
    else
        ks = Vector{Nothing}()
    end

    zs = [zero(u) for _ in 1:(s - 1)]
    push!(zs, nlsolver.z)
    tmp_cache = _esdirk_tmp_cache(alg, u, rate_prototype, uEltypeNoUnits)

    return ESDIRKIMEXCache(
        u, uprev, fsalfirst, zs, ks, tmp_cache, nlsolver, tab, alg.step_limiter!
    )
end

_nlsolver_γc(alg, tab) = (tab.Ai[tab.s, tab.s], tab.nlsolver_init_c)
_nlsolver_γc(::Trapezoid, tab) = (1 // 2, 1)
_nlsolver_γc(::ImplicitEuler, tab) = (1, 1)

function alg_cache(
        alg::_PureSDIRKAlg, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = ESDIRKIMEXTableau(alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = _nlsolver_γc(alg, tab)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    uprev3 = (alg isa Trapezoid) ? u : nothing
    tprev2 = (alg isa Trapezoid) ? t : nothing
    return ESDIRKIMEXConstantCache(nlsolver, tab, uprev3, tprev2)
end

function alg_cache(
        alg::_PureSDIRKAlg, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = ESDIRKIMEXTableau(alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = _nlsolver_γc(alg, tab)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)
    s = tab.s
    ks = f isa SplitFunction ? [zero(u) for _ in 1:s] : Vector{Nothing}()
    zs = [zero(u) for _ in 1:(s - 1)]
    push!(zs, nlsolver.z)
    tmp_cache = _esdirk_tmp_cache(alg, u, rate_prototype, uEltypeNoUnits)
    algebraic_vars = if (alg isa ImplicitEuler) && f.mass_matrix !== I
        # find_algebraic_vars_eqs is GPU-safe (broadcast-based, Diagonal-aware),
        # unlike `eachcol` which triggers scalar indexing on GPU arrays.
        find_algebraic_vars_eqs(f.mass_matrix)[1]
    else
        nothing
    end
    uprev3 = (alg isa Trapezoid) ? zero(u) : nothing
    tprev2 = (alg isa Trapezoid) ? t : nothing
    uprev2_field = alg_extrapolates(alg) ? uprev2 : nothing
    return ESDIRKIMEXCache(
        u, uprev, fsalfirst, zs, ks, tmp_cache, nlsolver, tab, _esdirk_step_limiter!(alg),
        uprev2_field, algebraic_vars, uprev3, tprev2
    )
end
