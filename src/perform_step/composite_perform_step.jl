function initialize!(integrator, cache::DefaultCache)
    cache.current = cache.choice_function(integrator)
    algs = integrator.alg.algs
    if cache.current == 1
        if !isdefined(cache, :cache1)
            cache.cache1 = alg_cache(algs[1], cache.args...)
        end
        initialize!(integrator, cache.cache1)
    elseif cache.current == 2
        if !isdefined(cache, :cache2)
            cache.cache2 = alg_cache(algs[2], cache.args...)
        end
        initialize!(integrator, cache.cache2)
        # the controller was initialized by default for algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, algs[1], algs[2])
    elseif cache.current == 3
        if !isdefined(cache, :cache3)
            cache.cache3 = alg_cache(algs[3], cache.args...)
        end
        initialize!(integrator, cache.cache3)
        # the controller was initialized by default for algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, algs[1], algs[3])
    elseif cache.current == 4
        if !isdefined(cache, :cache4)
            cache.cache4 = alg_cache(algs[4], cache.args...)
        end
        initialize!(integrator, cache.cache4)
        # the controller was initialized by default for algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, algs[1], algs[4])
    elseif cache.current == 5
        if !isdefined(cache, :cache5)
            cache.cache5 = alg_cache(algs[5], cache.args...)
        end
        initialize!(integrator, cache.cache5)
        # the controller was initialized by default for algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, algs[1], algs[5])
    elseif cache.current == 6
        if !isdefined(cache, :cache6)
            cache.cache6 = alg_cache(algs[6], cache.args...)
        end
        initialize!(integrator, cache.cache6)
        # the controller was initialized by default for algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, algs[1], algs[6])
    end
    resize!(integrator.k, integrator.kshortsize)
end

function initialize!(integrator, cache::CompositeCache)
    cache.current = cache.choice_function(integrator)
    if cache.current == 1
        initialize!(integrator, @inbounds(cache.caches[1]))
    elseif cache.current == 2
        initialize!(integrator, @inbounds(cache.caches[2]))
        # the controller was initialized by default for integrator.alg.algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, integrator.alg.algs[1],
            integrator.alg.algs[2])
    else
        initialize!(integrator, @inbounds(cache.caches[cache.current]))
        reset_alg_dependent_opts!(integrator.opts.controller, integrator.alg.algs[1],
            integrator.alg.algs[cache.current])
    end
    resize!(integrator.k, integrator.kshortsize)
end

function initialize!(integrator, cache::CompositeCache{Tuple{T1, T2}, F}) where {T1, T2, F}
    cache.current = cache.choice_function(integrator)
    if cache.current == 1
        initialize!(integrator, @inbounds(cache.caches[1]))
    elseif cache.current == 2
        initialize!(integrator, @inbounds(cache.caches[2]))
        reset_alg_dependent_opts!(integrator.opts.controller, integrator.alg.algs[1],
            integrator.alg.algs[2])
    end
    resize!(integrator.k, integrator.kshortsize)
end

"""
If the user mixes adaptive and non-adaptive algorithms then, right after
initialize!, make integrator.opts match the default adaptivity such that
the behaviour is consistent.
In particular, prevents dt ⟶ 0 if starting with non-adaptive alg and opts.adaptive=true,
and dt=cst if starting with adaptive alg and opts.adaptive=false.
"""
function ensure_behaving_adaptivity!(integrator, cache::Union{DefaultCache, CompositeCache})
    if anyadaptive(integrator.alg) && !isadaptive(integrator.alg)
        integrator.opts.adaptive = isadaptive(integrator.alg.algs[cache.current])
    end
end

function perform_step!(integrator, cache::DefaultCache, repeat_step = false)
    if cache.current == 1
        perform_step!(integrator, @inbounds(cache.cache1), repeat_step)
    elseif cache.current == 2
        perform_step!(integrator, @inbounds(cache.cache2), repeat_step)
    elseif cache.current == 3
        perform_step!(integrator, @inbounds(cache.cache3), repeat_step)
    elseif cache.current == 4
        perform_step!(integrator, @inbounds(cache.cache4), repeat_step)
    elseif cache.current == 5
        perform_step!(integrator, @inbounds(cache.cache5), repeat_step)
    elseif cache.current == 6
        perform_step!(integrator, @inbounds(cache.cache6), repeat_step)
    end
end

function perform_step!(integrator, cache::CompositeCache, repeat_step = false)
    if cache.current == 1
        perform_step!(integrator, @inbounds(cache.caches[1]), repeat_step)
    elseif cache.current == 2
        perform_step!(integrator, @inbounds(cache.caches[2]), repeat_step)
    else
        perform_step!(integrator, @inbounds(cache.caches[cache.current]), repeat_step)
    end
end

choose_algorithm!(integrator, cache::OrdinaryDiffEqCache) = nothing

function choose_algorithm!(integrator,
        cache::CompositeCache{Tuple{T1, T2}, F}) where {T1, T2, F}
    new_current = cache.choice_function(integrator)
    old_current = cache.current
    @inbounds if new_current != old_current
        cache.current = new_current
        if new_current == 1
            initialize!(integrator, @inbounds(cache.caches[1]))
        elseif new_current == 2
            initialize!(integrator, @inbounds(cache.caches[2]))
        end
        if old_current == 1 && new_current == 2
            reset_alg_dependent_opts!(integrator, integrator.alg.algs[1],
                integrator.alg.algs[2])
            transfer_cache!(integrator, integrator.cache.caches[1],
                integrator.cache.caches[2])
        elseif old_current == 2 && new_current == 1
            reset_alg_dependent_opts!(integrator, integrator.alg.algs[2],
                integrator.alg.algs[1])
            transfer_cache!(integrator, integrator.cache.caches[2],
                integrator.cache.caches[1])
        end
    end
end

function choose_algorithm!(integrator, cache::CompositeCache)
    new_current = cache.choice_function(integrator)
    old_current = cache.current
    @inbounds if new_current != old_current
        cache.current = new_current
        initialize!(integrator, @inbounds(cache.caches[new_current]))

        controller.beta2 = beta2_default(alg2)
        controller.beta1 = beta2_default(alg2)
        DEFAULTBETA2S

        reset_alg_dependent_opts!(integrator, integrator.alg.algs[old_current],
        integrator.alg.algs[new_current])
        transfer_cache!(integrator, integrator.cache.caches[old_current],
        integrator.cache.caches[new_current])
    end
end

function choose_algorithm!(integrator, cache::CompositeCache{<:Any, <:AutoSwitchCache{DefaultODESolver}})
    new_current = cache.choice_function(integrator)
    old_current = cache.current
    @inbounds if new_current != old_current
        cache.current = new_current
        if new_current == 1
            initialize!(integrator, @inbounds(cache.caches[1]))
        elseif new_current == 2
            initialize!(integrator, @inbounds(cache.caches[2]))
        elseif new_current == 3
            initialize!(integrator, @inbounds(cache.caches[3]))
        elseif new_current == 4
            initialize!(integrator, @inbounds(cache.caches[4]))
        elseif new_current == 5
            initialize!(integrator, @inbounds(cache.caches[5]))
        elseif new_current == 6
            initialize!(integrator, @inbounds(cache.caches[6]))
        else
            initialize!(integrator, @inbounds(cache.caches[new_current]))
        end

        # dtchangable, qmin_default, qmax_default, and isadaptive ignored since all same
        integrator.opts.controller.beta1 = DEFAULTBETA1S[new_current]
        integrator.opts.controller.beta2 = DEFAULTBETA2S[new_current]
    end
    nothing
end

"""
If no user default, then this will change the default to the defaults
for the second algorithm.
Except if the user default turns out to be the default for the current alg,
then it will change anyway and keep changing afterwards (e.g. adaptive).
"""
function reset_alg_dependent_opts!(integrator, alg1, alg2)
    integrator.dtchangeable = isdtchangeable(alg2)
    if integrator.opts.adaptive == isadaptive(alg1)
        integrator.opts.adaptive = isadaptive(alg2)
    end
    if integrator.opts.qmin == qmin_default(alg1)
        integrator.opts.qmin = qmin_default(alg2)
    end
    if integrator.opts.qmax == qmax_default(alg1)
        integrator.opts.qmax == qmax_default(alg2)
    end
    reset_alg_dependent_opts!(integrator.opts.controller, alg1, alg2)
    nothing
end

# Write how to transfer the cache variables from one cache to the other
# Example: send the history variables from one multistep method to another

transfer_cache!(integrator, alg1, alg2) = nothing
