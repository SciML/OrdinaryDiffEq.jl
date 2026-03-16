@inline function DiffEqNoiseProcess.setup_next_step!(integrator::SDEIntegrator)
    !isnothing(integrator.W) &&
        DiffEqNoiseProcess.setup_next_step!(integrator.W, integrator.u, integrator.p)
    return !isnothing(integrator.P) &&
        DiffEqNoiseProcess.setup_next_step!(integrator.P, integrator.u, integrator.p)
end

@inline function handle_callback_modifiers!(integrator::SDEIntegrator)
    #integrator.reeval_fsal = true
    return if integrator.P !== nothing && integrator.opts.adaptive
        if integrator.cache isa StochasticDiffEqMutableCache
            oldrate = integrator.P.cache.currate
            integrator.P.cache.rate(oldrate, integrator.u, integrator.p, integrator.t)
        else
            integrator.P.cache.currate = integrator.P.cache.rate(integrator.u, integrator.p, integrator.t)
        end
    end
end

@inline initialize!(integrator, cache::StochasticDiffEqCache, f = integrator.f) = nothing

function nlsolve!(integrator, cache)
    return DiffEqBase.nlsolve!(cache.nlsolver, cache.nlsolver.cache, integrator)
end

# TauLeapingDrift: wrapper for tau-leaping drift function used by nlsolver
# Computes drift(u, p, t) = c(u, p, t, rate(u, p, t), nothing)
# where c is the stoichiometry function and rate is the propensity function
struct TauLeapingDrift{C, R, RateCache, IIP}
    c::C              # Stoichiometry function (from integrator.c)
    rate::R           # Rate function (from integrator.P.cache.rate)
    rate_cache::RateCache  # Cache for rate values (for in-place version)
end

# Out-of-place: drift(u, p, t)
function (td::TauLeapingDrift{C, R, Nothing, false})(u, p, t) where {C, R}
    rates = td.rate(u, p, t)
    return td.c(u, p, t, rates, nothing)
end

# In-place: drift(du, u, p, t)
function (td::TauLeapingDrift{C, R, RateCache, true})(du, u, p, t) where {C, R, RateCache}
    td.rate(td.rate_cache, u, p, t)
    td.c(du, u, p, t, td.rate_cache, nothing)
    return nothing
end

# nlsolve_f override for ImplicitTauLeaping: return TauLeapingDrift wrapper
# Uses ODEIntegrator{A} directly since parameterized type aliases don't work.
function OrdinaryDiffEqCore.nlsolve_f(integrator::OrdinaryDiffEqCore.ODEIntegrator{A}) where {A <: ImplicitTauLeaping}
    # Determine if the cache is in-place or constant (out-of-place) based on cache type
    cache = integrator.cache
    if cache isa ImplicitTauLeapingCache
        # In-place version - use rate cache from the integrator's cache
        rate_cache = cache.rate_at_uprev
        return TauLeapingDrift{typeof(integrator.c), typeof(integrator.P.cache.rate), typeof(rate_cache), true}(
            integrator.c, integrator.P.cache.rate, rate_cache
        )
    else
        # Out-of-place (constant cache) version
        return TauLeapingDrift{typeof(integrator.c), typeof(integrator.P.cache.rate), Nothing, false}(
            integrator.c, integrator.P.cache.rate, nothing
        )
    end
end

# nlsolve_f override for ThetaTrapezoidalTauLeaping: return TauLeapingDrift wrapper
# Uses ODEIntegrator{A} directly since parameterized type aliases don't work.
function OrdinaryDiffEqCore.nlsolve_f(integrator::OrdinaryDiffEqCore.ODEIntegrator{A}) where {A <: ThetaTrapezoidalTauLeaping}
    # Determine if the cache is in-place or constant (out-of-place) based on cache type
    cache = integrator.cache
    if cache isa ThetaTrapezoidalTauLeapingCache
        # In-place version - use rate cache from the integrator's cache
        rate_cache = cache.rate_at_uprev
        return TauLeapingDrift{typeof(integrator.c), typeof(integrator.P.cache.rate), typeof(rate_cache), true}(
            integrator.c, integrator.P.cache.rate, rate_cache
        )
    else
        # Out-of-place (constant cache) version
        return TauLeapingDrift{typeof(integrator.c), typeof(integrator.P.cache.rate), Nothing, false}(
            integrator.c, integrator.P.cache.rate, nothing
        )
    end
end

# ============================================================================
# Traits for SDE composite algorithms/caches
# ============================================================================

# Trait: SDE composite algorithm types
function OrdinaryDiffEqCore.is_composite_algorithm(
        alg::Union{StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm},
    )
    return true
end

# Trait: SDE composite cache type
OrdinaryDiffEqCore.is_composite_cache(cache::StochasticCompositeCache) = true

# ============================================================================
# Noise interface methods for OrdinaryDiffEqCore's noise functions
# These implement the accept_noise!/reject_noise!/save_noise!/noise_curt/
# is_noise_saveable interface that ODE's unified loop functions call.
# ============================================================================

function OrdinaryDiffEqCore.accept_noise!(W::DiffEqNoiseProcess.AbstractNoiseProcess, dt, u, p, setup)
    return DiffEqNoiseProcess.accept_step!(W, dt, u, p, setup)
end

function OrdinaryDiffEqCore.reject_noise!(W::DiffEqNoiseProcess.AbstractNoiseProcess, dt, u, p)
    return DiffEqNoiseProcess.reject_step!(W, dt, u, p)
end

function OrdinaryDiffEqCore.save_noise!(W::DiffEqNoiseProcess.AbstractNoiseProcess)
    return DiffEqNoiseProcess.save_noise!(W)
end

OrdinaryDiffEqCore.noise_curt(W::DiffEqNoiseProcess.AbstractNoiseProcess) = W.curt

OrdinaryDiffEqCore.is_noise_saveable(W::NoiseProcess) = true
OrdinaryDiffEqCore.is_noise_saveable(W::DiffEqNoiseProcess.AbstractNoiseProcess) = false

# ============================================================================
# is_constant_cache for SDE cache types (needed by ODE's change_t_via_interpolation!)
# ============================================================================

OrdinaryDiffEqCore.is_constant_cache(::StochasticDiffEqConstantCache) = true
OrdinaryDiffEqCore.is_constant_cache(::StochasticDiffEqMutableCache) = false
OrdinaryDiffEqCore.is_constant_cache(cache::StochasticCompositeCache) =
    OrdinaryDiffEqCore.is_constant_cache(cache.caches[1])

# ============================================================================
# reinit_noise!: reinitialize noise process (called from ODE's reinit!)
# ============================================================================

function OrdinaryDiffEqCore.reinit_noise!(W::DiffEqNoiseProcess.AbstractNoiseProcess, dt)
    return DiffEqNoiseProcess.reinit!(W, dt)
end

# _determine_initdt: SDE extension (called from ODE's auto_dt_reset!)
function OrdinaryDiffEqCore._determine_initdt(integrator::SDEIntegrator)
    return OrdinaryDiffEqCore.sde_determine_initdt(
        integrator.u, integrator.t,
        integrator.tdir, integrator.opts.dtmax,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, integrator.sol.prob,
        get_current_alg_order(integrator.alg, integrator.cache),
        integrator
    )
end
