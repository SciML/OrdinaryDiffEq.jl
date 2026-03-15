@inline function initialize!(integrator, cache::StochasticCompositeCache)
    cache.current = cache.choice_function(integrator)
    return if cache.current == 1
        initialize!(integrator, @inbounds(cache.caches[1]))
    elseif cache.current == 2
        initialize!(integrator, @inbounds(cache.caches[2]))
    else
        initialize!(integrator, @inbounds(cache.caches[cache.current]))
    end
end

@inline function perform_step!(integrator, cache::StochasticCompositeCache)
    return if cache.current == 1
        perform_step!(integrator, @inbounds(cache.caches[1]))
    elseif cache.current == 2
        perform_step!(integrator, @inbounds(cache.caches[2]))
    else
        perform_step!(integrator, @inbounds(cache.caches[cache.current]))
    end
end

# choose_algorithm!, reset_alg_dependent_opts!, and transfer_cache! are now
# provided by OrdinaryDiffEqCore's generic composite infrastructure.
# The generic choose_algorithm!(integrator, cache) fallback in ODE uses the
# is_composite_cache trait to handle StochasticCompositeCache.
