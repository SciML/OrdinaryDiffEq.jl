abstract type SDIRKMutableCache <: OrdinaryDiffEqMutableCache end
abstract type SDIRKConstantCache <: OrdinaryDiffEqConstantCache end
function get_fsalfirstlast(cache::SDIRKMutableCache, u)
    return (cache.fsalfirst, du_alias_or_new(cache.nlsolver, cache.fsalfirst))
end

@cache mutable struct ImplicitEulerCache{
        uType, rateType, uNoUnitsType, N, AV, StepLimiter,
    } <:
    SDIRKMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    fsalfirst::rateType
    atmp::uNoUnitsType
    nlsolver::N
    algebraic_vars::AV
    step_limiter!::StepLimiter
end

function alg_cache(
        alg::ImplicitEuler, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    algebraic_vars = f.mass_matrix === I ? nothing :
        [all(iszero, x) for x in eachcol(f.mass_matrix)]

    return ImplicitEulerCache(
        u, uprev, uprev2, fsalfirst, atmp, nlsolver, algebraic_vars, alg.step_limiter!
    )
end

mutable struct ImplicitEulerConstantCache{N} <: SDIRKConstantCache
    nlsolver::N
end

function alg_cache(
        alg::ImplicitEuler, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    return ImplicitEulerConstantCache(nlsolver)
end

mutable struct TrapezoidConstantCache{uType, tType, N} <: SDIRKConstantCache
    uprev3::uType
    tprev2::tType
    nlsolver::N
end

function alg_cache(
        alg::Trapezoid, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1 // 2, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )

    uprev3 = u
    tprev2 = t

    return TrapezoidConstantCache(uprev3, tprev2, nlsolver)
end

@cache mutable struct TrapezoidCache{
        uType, rateType, uNoUnitsType, tType, N, StepLimiter,
    } <:
    SDIRKMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    fsalfirst::rateType
    atmp::uNoUnitsType
    uprev3::uType
    tprev2::tType
    nlsolver::N
    step_limiter!::StepLimiter
end

function alg_cache(
        alg::Trapezoid, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1 // 2, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    uprev3 = zero(u)
    tprev2 = t
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    return TrapezoidCache(
        u, uprev, uprev2, fsalfirst, atmp, uprev3, tprev2, nlsolver, alg.step_limiter!
    )
end

# Pure SDIRK types (non-IMEX) unified into ESDIRKIMEXCache
const _PureSDIRKAlg = Union{
    OrdinaryDiffEqNewtonNonAdaptiveSDIRKAlgorithm,
    OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm,
}

# step_limiter! accessor — only some pure SDIRK algorithms have the field
_esdirk_step_limiter!(alg::OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm) = alg.step_limiter!
_esdirk_step_limiter!(alg::Union{ImplicitMidpoint, SDIRK2, TRBDF2}) = alg.step_limiter!
_esdirk_step_limiter!(alg) = trivial_limiter!

# smooth_est accessor — only adaptive algorithms carry this flag
_esdirk_smooth_est(alg::OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm) = alg.smooth_est
_esdirk_smooth_est(alg::OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm) = alg.smooth_est
_esdirk_smooth_est(alg) = false

mutable struct ESDIRKIMEXConstantCache{Tab, N} <: OrdinaryDiffEqConstantCache
    nlsolver::N
    tab::Tab
end

mutable struct ESDIRKIMEXCache{uType, rateType, uNoUnitsType, N, Tab, kType, StepLimiter} <:
    SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    zs::Vector{uType}
    ks::Vector{kType}
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
    step_limiter!::StepLimiter
end

function full_cache(c::ESDIRKIMEXCache)
    base = (c.u, c.uprev, c.fsalfirst, c.zs..., c.atmp)
    if eltype(c.ks) !== Nothing
        return tuple(base..., c.ks...)
    end
    return base
end

function OrdinaryDiffEqCore.strip_cache(cache::ESDIRKIMEXCache)
    s = length(cache.zs)
    return SciMLBase.constructorof(typeof(cache))(
        nothing, nothing, nothing,
        Vector{Nothing}(undef, s),
        Vector{Nothing}(undef, s),
        nothing, nothing, nothing, nothing
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
        ks = Vector{Nothing}(nothing, s)
    end

    zs = [zero(u) for _ in 1:(s - 1)]
    push!(zs, nlsolver.z)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    return ESDIRKIMEXCache(
        u, uprev, fsalfirst, zs, ks, atmp, nlsolver, tab, alg.step_limiter!
    )
end

function alg_cache(
        alg::_PureSDIRKAlg, u, rate_prototype, ::Type{uEltypeNoUnits},
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
        alg::_PureSDIRKAlg, u, rate_prototype, ::Type{uEltypeNoUnits},
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
    ks = Vector{Nothing}(nothing, s)
    zs = [zero(u) for _ in 1:(s - 1)]
    push!(zs, nlsolver.z)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    return ESDIRKIMEXCache(
        u, uprev, fsalfirst, zs, ks, atmp, nlsolver, tab, _esdirk_step_limiter!(alg)
    )
end
