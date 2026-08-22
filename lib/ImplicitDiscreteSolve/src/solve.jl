# u === nothing path: nothing to step. The integrator's state is unchanged
# and the step trivially succeeds.
function perform_step!(
        integrator, cache::IDSolveCache{Nothing, Nothing}, repeat_step = false
    )
    return nothing
end

function perform_step!(integrator, cache::IDSolveCache, repeat_step = false)
    (; alg, u, uprev, dt, t, tprev, f, p) = integrator
    (; nlcache, observer) = cache

    # initial guess
    if alg.extrapolant == :constant
        cache.z .= integrator.u
    else
        error("Unknown extrapolant $(alg.extrapolant).")
    end
    state = ImplicitDiscreteState(cache.z, p, t + dt)

    # nonlinear solve step
    SciMLBase.reinit!(nlcache, cache.z; p = state)
    converged, znew = _solve_nonlinear!(nlcache, observer)
    if !converged
        integrator.force_stepfail = true
        return
    end

    # Accept step
    return u .= znew
end

"""
    _can_observe_steps(nlcache) -> Bool

Whether `nlcache` can be driven by `NonlinearSolveBase.solve_cache!` with a step observer.

Two cache types cannot. A `NonlinearSolveNoInitCache`, which is what the SimpleNonlinearSolve
algorithms return, holds the problem and options but no iteration state, so it has no `step!`
to observe; NonlinearSolve's own iterator-interface documentation gives `solve!` as the way to
drive it. A `NonlinearSolvePolyAlgorithmCache` does step, but `solve_cache!` finishes by
reading `cache.termination_cache` and `cache.trace`, which a polyalgorithm keeps on its active
subsolver rather than on itself (SciML/OrdinaryDiffEq.jl#4138). That second case is fixed
upstream by SciML/NonlinearSolve.jl#1189; this branch can go once that is released.
"""
@inline function _can_observe_steps(nlcache)
    return !(
        nlcache isa NonlinearSolveBase.NonlinearSolveNoInitCache ||
            nlcache isa NonlinearSolveBase.NonlinearSolvePolyAlgorithmCache
    )
end

# Returns `(converged, u)`. The state is returned rather than read back off the cache
# afterwards because only the observer path leaves the solution in `state_values`; a cache
# driven by `solve!` reports it on the returned solution and can leave its own state at the
# initial guess, which silently returns `u0` as the step (#4138).
function _solve_nonlinear!(nlcache, observer)
    reset!(observer)
    if _can_observe_steps(nlcache)
        retcode = NonlinearSolveBase.solve_cache!(nlcache; step_observer = observer)
        return retcode == ReturnCode.Success, state_values(nlcache)
    end
    # No step observer means no per-iteration convergence rates, so `Θks` stays empty.
    # `KantorovichTypeController` already copes: `stepsize_controller!` falls back to
    # `Θmin`, `step_reject_controller!` iterates nothing, and `accept_step_controller`
    # reduces to `all` over an empty collection, which is its documented `strict = false`
    # behaviour of accepting whenever the solve converged.
    sol = solve!(nlcache)
    return sol.retcode == ReturnCode.Success, sol.u
end

function initialize!(integrator, cache::IDSolveCache)
    return integrator.u isa AbstractVector && (cache.z .= integrator.u)
end

function _initialize_dae!(
        integrator, prob::ImplicitDiscreteProblem,
        alg::DefaultInit, x::Union{Val{true}, Val{false}}
    )
    isnothing(prob.u0) && return
    atol = one(eltype(prob.u0)) * 1.0e-12
    return if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(
            integrator, prob,
            OverrideInit(atol), x
        )
    else
        (; u, p, t) = integrator
        (; z, nlcache, observer) = integrator.cache
        z .= u
        initstate = ImplicitDiscreteState(z, p, t)
        SciMLBase.reinit!(nlcache, u; p = initstate)
        converged, unew = _solve_nonlinear!(nlcache, observer)
        if converged
            integrator.u .= unew
        else
            integrator.sol = SciMLBase.solution_new_retcode(
                integrator.sol, ReturnCode.InitialFailure
            )
        end
    end
end
