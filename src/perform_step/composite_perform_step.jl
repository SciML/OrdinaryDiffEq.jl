function initialize!(integrator, cache::CompositeCache{Fallbacks}) where Fallbacks
    cache.current = cache.choice_function(integrator)
    if cache.current == 1
        initialize!(integrator, @inbounds(cache.caches[1]))
    elseif cache.current == 2
        initialize!(integrator, @inbounds(cache.caches[2]))
        # the controller was initialized by default for integrator.alg.algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, integrator.alg.algs[1],
            integrator.alg.algs[2])
    elseif cache.current == 3
        initialize!(integrator, @inbounds(cache.caches[3]))
        # the controller was initialized by default for integrator.alg.algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, integrator.alg.algs[1],
            integrator.alg.algs[3])
    elseif cache.current == 4
        initialize!(integrator, @inbounds(cache.caches[4]))
        # the controller was initialized by default for integrator.alg.algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, integrator.alg.algs[1],
            integrator.alg.algs[4])
    elseif cache.current == 5
        initialize!(integrator, @inbounds(cache.caches[5]))
        # the controller was initialized by default for integrator.alg.algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, integrator.alg.algs[1],
            integrator.alg.algs[5])
    elseif cache.current == 6
        initialize!(integrator, @inbounds(cache.caches[6]))
        # the controller was initialized by default for integrator.alg.algs[1]
        reset_alg_dependent_opts!(integrator.opts.controller, integrator.alg.algs[1],
            integrator.alg.algs[6])
    elseif Fallbacks
        initialize!(integrator, @inbounds(cache.caches[cache.current]))
        reset_alg_dependent_opts!(integrator.opts.controller, integrator.alg.algs[1],
            integrator.alg.algs[cache.current])
    end
    resize!(integrator.k, integrator.kshortsize)
end

function initialize!(integrator, cache::CompositeCache{false, Tuple{T1, T2}, F}) where {T1, T2, F}
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
function ensure_behaving_adaptivity!(integrator, cache::CompositeCache)
    if anyadaptive(integrator.alg) && !isadaptive(integrator.alg)
        integrator.opts.adaptive = isadaptive(integrator.alg.algs[cache.current])
    end
end

function perform_step!(integrator, cache::CompositeCache{Fallbacks}, repeat_step = false) where Fallbacks
    if cache.current == 1
        perform_step!(integrator, @inbounds(cache.caches[1]), repeat_step)
    elseif cache.current == 2
        perform_step!(integrator, @inbounds(cache.caches[2]), repeat_step)
    elseif cache.current == 3
        perform_step!(integrator, @inbounds(cache.caches[3]), repeat_step)
    elseif cache.current == 4
        perform_step!(integrator, @inbounds(cache.caches[4]), repeat_step)
    elseif cache.current == 5
        perform_step!(integrator, @inbounds(cache.caches[5]), repeat_step)
    elseif cache.current == 6
        perform_step!(integrator, @inbounds(cache.caches[6]), repeat_step)
    elseif Fallbacks
        perform_step!(integrator, @inbounds(cache.caches[cache.current]), repeat_step)
    end
end

function perform_step!(integrator, cache::CompositeCache{false, Tuple{T1, T2}, F},
    repeat_step = false) where {T1, T2, F}
    if cache.current == 1
        perform_step!(integrator, @inbounds(cache.caches[1]), repeat_step)
    elseif cache.current == 2
        perform_step!(integrator, @inbounds(cache.caches[2]), repeat_step)
    end
end

choose_algorithm!(integrator, cache::OrdinaryDiffEqCache) = nothing

function choose_algorithm!(integrator,
    cache::CompositeCache{false, Tuple{T1, T2}, F}) where {T1, T2, F}
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

function choose_algorithm!(integrator, cache::CompositeCache{Fallbacks, T, F}) where {Fallbacks, T, F}
    new_current = cache.choice_function(integrator)
    old_current = cache.current
    !Fallbacks && error("Hitting fallbacks!")
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

function choose_algorithm!(integrator, cache::CompositeCache{Fallbacks, T, <:AutoSwitchCache{DefaultODESolver}}) where {Fallbacks, T}
    new_current = cache.choice_function(integrator)
    old_current = cache.current
    @inbounds if new_current != old_current
        cache.current = new_current
        if new_current == 1
            initialize!(integrator, @inbounds(cache.caches[1])); nothing
        elseif new_current == 2
            initialize!(integrator, @inbounds(cache.caches[2])); nothing
        elseif new_current == 3
            initialize!(integrator, @inbounds(cache.caches[3])); nothing
        elseif new_current == 4
            initialize!(integrator, @inbounds(cache.caches[4])); nothing
        elseif new_current == 5
            initialize!(integrator, @inbounds(cache.caches[5])); nothing
        elseif new_current == 6
            initialize!(integrator, @inbounds(cache.caches[6])); nothing
        else
            error("Unrachable reached. Report this error")
        end

        # dtchangable, qmin_default, qmax_default, and isadaptive ignored since all same
        integrator.opts.controller.beta1 = DEFAULTBETA1S[new_current]
        integrator.opts.controller.beta2 = DEFAULTBETA2S[new_current]
    end
end

#=
"""
function choose_algorithm!(integrator, cache::CompositeCache{Fallbacks}) where Fallbacks
    new_current = cache.choice_function(integrator)
    old_current = cache.current
    @inbounds if new_current != old_current
        cache.current = new_current
        initialize!(integrator, @inbounds(cache.caches[new_current]))
        reset_alg_dependent_opts!(integrator, integrator.alg.algs[old_current],
        integrator.alg.algs[new_current])
        transfer_cache!(integrator, integrator.cache.caches[old_current],
        integrator.cache.caches[new_current])
    end
end
"""
@generated function choose_algorithm!(integrator, cache::CompositeCache{Fallbacks, T, F}) where {Fallbacks, T, F}
    initialize_ex = :()
    for idx in 1:length(T.types)
        newex = quote
            initialize!(integrator, @inbounds(cache.caches[$idx]))
        end
        initialize_ex = if initialize_ex == :()
            Expr(:elseif, :(new_current == $idx)), newex,
                :(error("Algorithm Choice not Allowed"))
        else
            Expr(:elseif, :(new_current == $idx), newex, initialize_ex)
        end
    end
    initialize_ex = Expr(:if, initialize_ex.args...)

    swap_ex = :()
    for idx in 1:length(T.types), idx2 in 1:length(T.types)
        new_swap_ex = quote
            reset_alg_dependent_opts!(integrator, integrator.alg.algs[$idx],
            integrator.alg.algs[$idx2])
            transfer_cache!(integrator, integrator.cache.caches[$idx],
            integrator.cache.caches[$idx2])
        end
        swap_ex = if swap_ex == :()
            Expr(:elseif, :(old_current == $idx && new_current == $idx2)), new_swap_ex,
                :(error("Algorithm Choice not Allowed"))
        else
            Expr(:elseif, :(new_current == $idx), swap_ex, swap_ex)
        end
    end
    swap_ex = Expr(:if, swap_ex.args...)


    quote
        new_current = cache.choice_function(integrator)
        old_current = cache.current
        @inbounds if new_current != old_current
            cache.current = new_current
            $initialize_ex
            $swap_ex
        end
    end
end
=#

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
