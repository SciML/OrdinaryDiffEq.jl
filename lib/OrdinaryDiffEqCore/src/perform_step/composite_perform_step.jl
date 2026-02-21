function init_ith_default_cache(cache::DefaultCacheType, algs, i)
    return if i == 1
        if !isdefined(cache, :cache1)
            cache.cache1 = alg_cache(algs[1], cache.args...)
        end
    elseif i == 2
        if !isdefined(cache, :cache2)
            cache.cache2 = alg_cache(algs[2], cache.args...)
        end
    elseif i == 3
        if !isdefined(cache, :cache3)
            cache.cache3 = alg_cache(algs[3], cache.args...)
        end
    elseif i == 4
        if !isdefined(cache, :cache4)
            cache.cache4 = alg_cache(algs[4], cache.args...)
        end
    elseif i == 5
        if !isdefined(cache, :cache5)
            cache.cache5 = alg_cache(algs[5], cache.args...)
        end
    elseif i == 6
        if !isdefined(cache, :cache6)
            cache.cache6 = alg_cache(algs[6], cache.args...)
        end
    end
end

function initialize!(integrator, cache::DefaultCacheType)
    cache.current = cache.choice_function(integrator)
    algs = integrator.alg.algs
    init_ith_default_cache(cache, algs, cache.current)
    u = integrator.u
    if cache.current == 1
        fsalfirst, fsallast = get_fsalfirstlast(cache.cache1, u)
        !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
        !isnothing(fsallast) && (integrator.fsallast = fsallast)
        initialize!(integrator, cache.cache1)
    elseif cache.current == 2
        fsalfirst, fsallast = get_fsalfirstlast(cache.cache2, u)
        !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
        !isnothing(fsallast) && (integrator.fsallast = fsallast)
        initialize!(integrator, cache.cache2)
        # the controller was initialized by default for algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, algs[1], algs[2])
    elseif cache.current == 3
        fsalfirst, fsallast = get_fsalfirstlast(cache.cache3, u)
        !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
        !isnothing(fsallast) && (integrator.fsallast = fsallast)
        initialize!(integrator, cache.cache3)
        # the controller was initialized by default for algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, algs[1], algs[3])
    elseif cache.current == 4
        fsalfirst, fsallast = get_fsalfirstlast(cache.cache4, u)
        !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
        !isnothing(fsallast) && (integrator.fsallast = fsallast)
        initialize!(integrator, cache.cache4)
        # the controller was initialized by default for algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, algs[1], algs[4])
    elseif cache.current == 5
        fsalfirst, fsallast = get_fsalfirstlast(cache.cache5, u)
        !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
        !isnothing(fsallast) && (integrator.fsallast = fsallast)
        initialize!(integrator, cache.cache5)
        # the controller was initialized by default for algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, algs[1], algs[5])
    elseif cache.current == 6
        fsalfirst, fsallast = get_fsalfirstlast(cache.cache6, u)
        !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
        !isnothing(fsallast) && (integrator.fsallast = fsallast)
        initialize!(integrator, cache.cache6)
        # the controller was initialized by default for algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, algs[1], algs[6])
    end
    return resize!(integrator.k, integrator.kshortsize)
end

function initialize!(integrator, cache::CompositeCache)
    cache.current = cache.choice_function(integrator)
    u = integrator.u
    if cache.current == 1
        fsalfirst, fsallast = get_fsalfirstlast(cache.caches[1], u)
        !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
        !isnothing(fsallast) && (integrator.fsallast = fsallast)
        initialize!(integrator, @inbounds(cache.caches[1]))
    elseif cache.current == 2
        fsalfirst, fsallast = get_fsalfirstlast(cache.caches[2], u)
        !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
        !isnothing(fsallast) && (integrator.fsallast = fsallast)
        initialize!(integrator, @inbounds(cache.caches[2]))
        # the controller was initialized by default for integrator.alg.algs[1]
        reset_alg_dependent_opts!(
            integrator.opts.controller, integrator.alg.algs[1],
            integrator.alg.algs[2]
        )
    else
        fsalfirst, fsallast = get_fsalfirstlast(cache.caches[cache.current], u)
        !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
        !isnothing(fsallast) && (integrator.fsallast = fsallast)
        initialize!(integrator, @inbounds(cache.caches[cache.current]))
        reset_alg_dependent_opts!(
            integrator.opts.controller, integrator.alg.algs[1],
            integrator.alg.algs[cache.current]
        )
    end
    return resize!(integrator.k, integrator.kshortsize)
end

function initialize!(integrator, cache::CompositeCache{Tuple{T1, T2}, F}) where {T1, T2, F}
    cache.current = cache.choice_function(integrator)
    u = integrator.u
    if cache.current == 1
        fsalfirst, fsallast = get_fsalfirstlast(cache.caches[1], u)
        !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
        !isnothing(fsallast) && (integrator.fsallast = fsallast)
        initialize!(integrator, @inbounds(cache.caches[1]))
    elseif cache.current == 2
        fsalfirst, fsallast = get_fsalfirstlast(cache.caches[2], u)
        !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
        !isnothing(fsallast) && (integrator.fsallast = fsallast)
        initialize!(integrator, @inbounds(cache.caches[2]))
        reset_alg_dependent_opts!(
            integrator.opts.controller, integrator.alg.algs[1],
            integrator.alg.algs[2]
        )
    end
    return resize!(integrator.k, integrator.kshortsize)
end

"""
If the user mixes adaptive and non-adaptive algorithms then, right after
initialize!, make integrator.opts match the default adaptivity such that
the behaviour is consistent.
In particular, prevents dt ⟶ 0 if starting with non-adaptive alg and opts.adaptive=true,
and dt=cst if starting with adaptive alg and opts.adaptive=false.
"""
function ensure_behaving_adaptivity!(integrator, cache::Union{DefaultCacheType, CompositeCache})
    return if anyadaptive(integrator.alg) && !isadaptive(integrator.alg)
        integrator.opts.adaptive = isadaptive(integrator.alg.algs[cache.current])
    end
end

function perform_step!(integrator, cache::DefaultCacheType, repeat_step = false)
    algs = integrator.alg.algs
    init_ith_default_cache(cache, algs, cache.current)
    return if cache.current == 1
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
    return if cache.current == 1
        perform_step!(integrator, @inbounds(cache.caches[1]), repeat_step)
    elseif cache.current == 2
        perform_step!(integrator, @inbounds(cache.caches[2]), repeat_step)
    else
        perform_step!(integrator, @inbounds(cache.caches[cache.current]), repeat_step)
    end
end

choose_algorithm!(integrator, cache::OrdinaryDiffEqCache) = nothing

function choose_algorithm!(
        integrator,
        cache::CompositeCache{Tuple{T1, T2}, F}
    ) where {T1, T2, F}
    new_current = cache.choice_function(integrator)
    old_current = cache.current
    u = integrator.u
    return @inbounds if new_current != old_current
        cache.current = new_current
        if new_current == 1
            fsalfirst, fsallast = get_fsalfirstlast(cache.caches[1], u)
            !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
            !isnothing(fsallast) && (integrator.fsallast = fsallast)
            initialize!(integrator, @inbounds(cache.caches[1]))
        elseif new_current == 2
            fsalfirst, fsallast = get_fsalfirstlast(cache.caches[2], u)
            !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
            !isnothing(fsallast) && (integrator.fsallast = fsallast)
            initialize!(integrator, @inbounds(cache.caches[2]))
        end
        if old_current == 1 && new_current == 2
            reset_alg_dependent_opts!(
                integrator, integrator.alg.algs[1],
                integrator.alg.algs[2]
            )
            transfer_cache!(
                integrator, integrator.cache.caches[1],
                integrator.cache.caches[2]
            )
        elseif old_current == 2 && new_current == 1
            reset_alg_dependent_opts!(
                integrator, integrator.alg.algs[2],
                integrator.alg.algs[1]
            )
            transfer_cache!(
                integrator, integrator.cache.caches[2],
                integrator.cache.caches[1]
            )
        end
    end
end

function choose_algorithm!(integrator, cache::DefaultCacheType)
    new_current = cache.choice_function(integrator)
    old_current = cache.current
    u = integrator.u
    return @inbounds if new_current != old_current
        algs = integrator.alg.algs
        cache.current = new_current
        init_ith_default_cache(cache, algs, new_current)
        new_cache = nothing
        old_cache = nothing
        if new_current == 1
            fsalfirst, fsallast = get_fsalfirstlast(cache.cache1, u)
            !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
            !isnothing(fsallast) && (integrator.fsallast = fsallast)
            initialize!(integrator, @inbounds(cache.cache1))
            new_cache = cache.cache1
        elseif new_current == 2
            fsalfirst, fsallast = get_fsalfirstlast(cache.cache2, u)
            !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
            !isnothing(fsallast) && (integrator.fsallast = fsallast)
            initialize!(integrator, @inbounds(cache.cache2))
            new_cache = cache.cache2
        elseif new_current == 3
            fsalfirst, fsallast = get_fsalfirstlast(cache.cache3, u)
            !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
            !isnothing(fsallast) && (integrator.fsallast = fsallast)
            initialize!(integrator, @inbounds(cache.cache3))
            new_cache = cache.cache3
        elseif new_current == 4
            fsalfirst, fsallast = get_fsalfirstlast(cache.cache4, u)
            !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
            !isnothing(fsallast) && (integrator.fsallast = fsallast)
            initialize!(integrator, @inbounds(cache.cache4))
            new_cache = cache.cache4
        elseif new_current == 5
            fsalfirst, fsallast = get_fsalfirstlast(cache.cache5, u)
            !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
            !isnothing(fsallast) && (integrator.fsallast = fsallast)
            initialize!(integrator, @inbounds(cache.cache5))
            new_cache = cache.cache5
        elseif new_current == 6
            fsalfirst, fsallast = get_fsalfirstlast(cache.cache6, u)
            !isnothing(fsalfirst) && (integrator.fsalfirst = fsalfirst)
            !isnothing(fsallast) && (integrator.fsallast = fsallast)
            initialize!(integrator, @inbounds(cache.cache6))
            new_cache = cache.cache6
        end

        if old_current == 1
            old_cache = cache.cache1
        elseif old_current == 2
            old_cache = cache.cache2
        elseif old_current == 3
            old_cache = cache.cache3
        elseif old_current == 4
            old_cache = cache.cache4
        elseif old_current == 5
            old_cache = cache.cache5
        elseif old_current == 6
            old_cache = cache.cache6
        end

        reset_alg_dependent_opts!(integrator, algs[old_current], algs[new_current])
        transfer_cache!(integrator, old_cache, new_cache)
    end
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
    reset_alg_dependent_opts!(integrator.opts.controller, alg1, alg2)
    return nothing
end

# Write how to transfer the cache variables from one cache to the other
# Example: send the history variables from one multistep method to another

transfer_cache!(integrator, alg1, alg2) = nothing
