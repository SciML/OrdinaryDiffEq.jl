function perform_step!(integrator, cache::IDSolveCache, repeat_step = false)
    (; alg, u, uprev, dt, t, tprev, f, p) = integrator
    (; nlcache, Θks) = cache

    # initial guess
    if alg.extrapolant == :constant
        cache.z .= integrator.u
    else
        error("Unknown extrapolant $(alg.extrapolant).")
    end
    state = ImplicitDiscreteState(cache.z, p, t + dt)

    # nonlinear solve step
    SciMLBase.reinit!(nlcache, p = state)

    # solve!(nlcache)
    # The solve here is simply unrolled by hand to query the convergence rate estimates "manually" for now
    if nlcache.retcode == ReturnCode.InitialFailure
        integrator.force_stepfail = true
        return
    end

    resize!(Θks, 0)
    residualnormprev = zero(eltype(u))
    while NonlinearSolveBase.not_terminated(nlcache)
        step!(nlcache)
        residualnorm = NonlinearSolveBase.L2_NORM(nlcache.fu)
        if nlcache.nsteps > 1
            # Θk = min(residualnorm/residualnormprev, incrementnorm/incrementnormprev)
            Θk = residualnorm / residualnormprev
            if residualnormprev ≈ 0.0 #|| incrementnormprev ≈ 0.0
                push!(Θks, 0.0)
            else
                push!(Θks, Θk)
            end
            # if nlcache.parameters.enforce_monotonic_convergence && Θk ≥ 1.0
            #     @debug "Newton-Raphson diverged. Aborting. ||r|| = $residualnorm" _group=:nlsolve
            #     return false
            # end
        end
        residualnormprev = residualnorm
    end

    # The solver might have set a different `retcode`
    if nlcache.retcode == ReturnCode.Default
        nlcache.retcode = ifelse(
            nlcache.nsteps ≥ nlcache.maxiters, ReturnCode.MaxIters, ReturnCode.Success
        )
    end

    NonlinearSolveBase.update_from_termination_cache!(nlcache.termination_cache, nlcache)

    NonlinearSolveBase.update_trace!(
        nlcache.trace, nlcache.nsteps, NonlinearSolveBase.get_u(nlcache),
        NonlinearSolveBase.get_fu(nlcache), nothing, nothing, nothing;
        last = Val(true)
    )

    if nlcache.retcode != ReturnCode.Success
        integrator.force_stepfail = true
        return
    end

    # Accept step
    return u .= nlcache.u
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
        (; u, p, t, f) = integrator
        initstate = ImplicitDiscreteState(u, p, t)

        _f = if isinplace(f)
            (resid, u_next, p) -> f(resid, u_next, p.u, p.p, p.t)
        else
            (u_next, p) -> f(u_next, p.u, p.p, p.t)
        end

        nlls = !isnothing(f.resid_prototype) &&
            (length(f.resid_prototype) != length(integrator.u))
        prob = if nlls
            NonlinearLeastSquaresProblem{isinplace(f)}(
                NonlinearFunction(_f; resid_prototype = f.resid_prototype), u, initstate
            )
        else
            NonlinearProblem{isinplace(f)}(_f, u, initstate)
        end
        sol = solve(prob, integrator.alg.nlsolve)
        if sol.retcode == ReturnCode.Success
            integrator.u = sol
        else
            integrator.sol = SciMLBase.solution_new_retcode(
                integrator.sol, ReturnCode.InitialFailure
            )
        end
    end
end
