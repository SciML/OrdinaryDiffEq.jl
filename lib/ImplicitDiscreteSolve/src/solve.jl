"""
    nlsolve_tolerances(integrator) -> NamedTuple

The `abstol`/`reltol` the step's nonlinear solve should be run at, taken from the tolerances
the integrator was initialized with.

An `ImplicitDiscreteProblem` carries no tolerances of its own, so `OrdinaryDiffEqCore`
stores `false` for one that was not supplied; that maps to `nothing`, which restores the
nonlinear solver's own default. A per-component tolerance array is reduced to its tightest
entry, because the nonlinear solve takes a scalar.
"""
function nlsolve_tolerances(integrator)
    return (;
        abstol = nlsolve_tolerance(integrator.opts.abstol),
        reltol = nlsolve_tolerance(integrator.opts.reltol),
    )
end

nlsolve_tolerance(::Nothing) = nothing
nlsolve_tolerance(::Bool) = nothing
nlsolve_tolerance(tol::Number) = tol
nlsolve_tolerance(tol::AbstractArray) = isempty(tol) ? nothing : minimum(tol)

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
    SciMLBase.reinit!(nlcache, cache.z; p = state, nlsolve_tolerances(integrator)...)
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

Whether the nonlinear solve is driven step by step through `NonlinearSolveBase.solve_cache!`
with the step observer, or run to completion with `solve!`.

`NonlinearSolveNoInitCache` is NonlinearSolveBase's public marker for a solver without the
iterator interface, so it takes the `solve!` path. A polyalgorithm cache does too, even
though it can be stepped: observing it feeds every subsolver's iterations, including the
ones that fail before the next subsolver takes over, to the step controller, which drives
`dt` to zero on the first step.
"""
@inline _can_observe_steps(nlcache) = true
@inline _can_observe_steps(::NonlinearSolveBase.NonlinearSolveNoInitCache) = false
@inline _can_observe_steps(::NonlinearSolveBase.NonlinearSolvePolyAlgorithmCache) = false

function _solve_nonlinear!(nlcache, observer)
    reset!(observer)
    if _can_observe_steps(nlcache)
        retcode = NonlinearSolveBase.solve_cache!(nlcache; step_observer = observer)
        return retcode == ReturnCode.Success, state_values(nlcache)
    end
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
        SciMLBase.reinit!(nlcache, u; p = initstate, nlsolve_tolerances(integrator)...)
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
