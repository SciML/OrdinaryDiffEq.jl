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
    if !_solve_nonlinear!(nlcache, observer)
        integrator.force_stepfail = true
        return
    end

    # Accept step
    return u .= state_values(nlcache)
end

function _solve_nonlinear!(nlcache, observer)
    reset!(observer)
    retcode = NonlinearSolveBase.solve_cache!(nlcache; step_observer = observer)
    return retcode == ReturnCode.Success
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
        if _solve_nonlinear!(nlcache, observer)
            integrator.u .= state_values(nlcache)
        else
            integrator.sol = SciMLBase.solution_new_retcode(
                integrator.sol, ReturnCode.InitialFailure
            )
        end
    end
end
