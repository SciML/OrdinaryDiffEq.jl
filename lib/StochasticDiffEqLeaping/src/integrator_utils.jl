import OrdinaryDiffEqCore: nlsolve_f

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
