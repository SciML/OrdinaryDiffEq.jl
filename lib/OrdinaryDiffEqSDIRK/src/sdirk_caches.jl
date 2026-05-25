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

mutable struct ESDIRKIMEXConstantCache{Tab, N, U3, T2} <: OrdinaryDiffEqConstantCache
    nlsolver::N
    tab::Tab
    uprev3::U3
    tprev2::T2
end

function ESDIRKIMEXConstantCache(nlsolver, tab)
    return ESDIRKIMEXConstantCache(nlsolver, tab, nothing, nothing)
end

mutable struct ESDIRKIMEXCache{
        uType, rateType, uNoUnitsType, N, Tab, kType, StepLimiter, U2, AV, U3, T2,
    } <: SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    zs::Vector{uType}
    ks::Vector{kType}
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
    step_limiter!::StepLimiter
    uprev2::U2
    algebraic_vars::AV
    uprev3::U3
    tprev2::T2
end

function ESDIRKIMEXCache(u, uprev, fsalfirst, zs, ks, atmp, nlsolver, tab, step_limiter!)
    return ESDIRKIMEXCache(
        u, uprev, fsalfirst, zs, ks, atmp, nlsolver, tab, step_limiter!,
        nothing, nothing, nothing, nothing
    )
end

function full_cache(c::ESDIRKIMEXCache)
    base = (c.u, c.uprev, c.fsalfirst, c.zs..., c.atmp)
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
    return SciMLBase.constructorof(typeof(cache))(
        nothing, nothing, nothing,
        Vector{Nothing}(undef, s),
        Vector{Nothing}(undef, s),
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing
    )
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
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    return ESDIRKIMEXCache(
        u, uprev, fsalfirst, zs, ks, atmp, nlsolver, tab, alg.step_limiter!
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
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    algebraic_vars = if (alg isa ImplicitEuler) && f.mass_matrix !== I
        [all(iszero, x) for x in eachcol(f.mass_matrix)]
    else
        nothing
    end
    uprev3 = (alg isa Trapezoid) ? zero(u) : nothing
    tprev2 = (alg isa Trapezoid) ? t : nothing
    uprev2_field = alg_extrapolates(alg) ? uprev2 : nothing
    return ESDIRKIMEXCache(
        u, uprev, fsalfirst, zs, ks, atmp, nlsolver, tab, _esdirk_step_limiter!(alg),
        uprev2_field, algebraic_vars, uprev3, tprev2
    )
end
