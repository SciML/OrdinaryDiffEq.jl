function perform_step!(integrator, cache::IDSolveCache, repeat_step = false)
    (; alg, u, uprev, dt, t, f, p) = integrator

    # initial guess
    if alg.extrapolant == :linear
        @.. broadcast=false cache.z=integrator.uprev + dt * (integrator.uprev - integrator.uprev2)
    else # :constant
        cache.z .= integrator.u
    end
    state = ImplicitDiscreteState(cache.z, p, t)

    # nonlinear solve step
    SciMLBase.reinit!(cache.nlcache, p=state)
    # TODO compute convergence rate estimate
    # for i in 1:10
    #     step!(cache.nlcache)
    #     # ...
    # end
    solve!(cache.nlcache)
    if cache.nlcache.retcode != ReturnCode.Success
        integrator.force_stepfail = true
        return
    end

    # Accept step
    u .= cache.nlcache.u
end

function initialize!(integrator, cache::IDSolveCache)
    integrator.u isa AbstractVector && (cache.z .= integrator.u)
end

function _initialize_dae!(integrator, prob::ImplicitDiscreteProblem,
        alg::DefaultInit, x::Union{Val{true}, Val{false}})
    isnothing(prob.u0) && return
    atol = one(eltype(prob.u0)) * 1e-12
    if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(integrator, prob,
            OverrideInit(atol), x)
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
                NonlinearFunction(_f; resid_prototype = f.resid_prototype), u, initstate)
        else
            NonlinearProblem{isinplace(f)}(_f, u, initstate)
        end
        sol = solve(prob, integrator.alg.nlsolve)
        if sol.retcode == ReturnCode.Success
            integrator.u = sol
        else
            integrator.sol = SciMLBase.solution_new_retcode(integrator.sol, ReturnCode.InitialFailure)
        end
    end
end
