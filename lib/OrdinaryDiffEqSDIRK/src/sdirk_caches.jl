abstract type SDIRKMutableCache <: OrdinaryDiffEqMutableCache end
abstract type SDIRKConstantCache <: OrdinaryDiffEqConstantCache end
function get_fsalfirstlast(cache::SDIRKMutableCache, u)
    (cache.fsalfirst, du_alias_or_new(cache.nlsolver, cache.fsalfirst))
end

# Unified SDIRK caches that work with any SDIRK tableau
@cache mutable struct SDIRKCache{uType, rateType, uNoUnitsType, Tab, N, AV, StepLimiter} <: SDIRKMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    fsalfirst::rateType
    zs::Vector{uType}
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
    algebraic_vars::AV
    step_limiter!::StepLimiter
    # For algorithms that need additional history
    uprev3::uType
    tprev2::typeof(1.0)
    zprev::uType
    zᵧ::uType
end

mutable struct SDIRKConstantCache{N, Tab, uType, tType} <: OrdinaryDiffEqConstantCache
    nlsolver::N
    tab::Tab
    # For algorithms that need additional history
    uprev3::uType
    tprev2::tType
end

# Unified alg_cache functions for mutable cache
function alg_cache(alg::Union{ImplicitEuler, ImplicitMidpoint, Trapezoid, TRBDF2, SDIRK2, SDIRK22, SSPSDIRK2,
                            Cash4, SFSDIRK4, SFSDIRK5, SFSDIRK6, SFSDIRK7, SFSDIRK8,
                            Hairer4, Hairer42, ESDIRK54I8L2SA, ESDIRK436L2SA2, ESDIRK437L2SA,
                            ESDIRK547L2SA2, ESDIRK659L2SA,
                            Kvaerno3, KenCarp3, CFNLIRK3, Kvaerno4, Kvaerno5, KenCarp4, 
                            KenCarp47, KenCarp5, KenCarp58},
                   u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}

    tab = get_sdirk_tableau(alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.γ
    s = length(tab.b)

    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true))
    fsalfirst = zero(rate_prototype)

    # Initialize all z stage vectors
    zs = Vector{typeof(u)}(undef, s)
    for i in 1:s
        zs[i] = zero(u)
    end
    zs[end] = nlsolver.z # use nlsolver.z for the last stage

    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    # Handle mass matrix for algebraic variables (ImplicitEuler needs this)
    algebraic_vars = f.mass_matrix === I ? nothing :
                     [all(iszero, x) for x in eachcol(f.mass_matrix)]

    # Get step limiter
    step_limiter! = hasproperty(alg, :step_limiter!) ? alg.step_limiter! : trivial_limiter!

    # Additional variables for algorithms that need history
    uprev3 = zero(u)
    tprev2 = t
    zprev = zero(u)
    zᵧ = zero(u)

    SDIRKCache(u, uprev, uprev2, fsalfirst, zs,
               atmp, nlsolver, tab, algebraic_vars, step_limiter!, uprev3, tprev2, zprev, zᵧ)
end

# Unified alg_cache functions for constant cache
function alg_cache(alg::Union{ImplicitEuler, ImplicitMidpoint, Trapezoid, TRBDF2, SDIRK2, SDIRK22, SSPSDIRK2,
                            Cash4, SFSDIRK4, SFSDIRK5, SFSDIRK6, SFSDIRK7, SFSDIRK8,
                            Hairer4, Hairer42, ESDIRK54I8L2SA, ESDIRK436L2SA2, ESDIRK437L2SA,
                            ESDIRK547L2SA2, ESDIRK659L2SA,
                            Kvaerno3, KenCarp3, CFNLIRK3, Kvaerno4, Kvaerno5, KenCarp4, 
                            KenCarp47, KenCarp5, KenCarp58},
                   u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                   uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}

    tab = get_sdirk_tableau(alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.γ

    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false))

    # Additional variables for algorithms that need history
    uprev3 = u
    tprev2 = t

    SDIRKConstantCache(nlsolver, tab, uprev3, tprev2)
end

# Keep old caches for backward compatibility for now, will be removed later.
