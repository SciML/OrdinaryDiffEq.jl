# Remake the nonlinear problem, then update
function perform_step!(integrator, cache::IDSolveCache, repeat_step = false)
    (; alg, u, uprev, dt, t, f, p) = integrator
    (; state, prob) = cache
    state.u .= uprev
    state.t_next = t
    prob = remake(prob, p = state)

    u = solve(prob, SimpleNewtonRaphson())
    integrator.sol = SciMLBase.solution_new_retcode(integrator.sol, u.retcode)
    integrator.u = u
end

function initialize!(integrator, cache::IDSolveCache)
    integrator.u isa AbstractVector && (cache.state.u .= integrator.u)
    cache.state.p = integrator.p
    cache.state.t_next = integrator.t
    f = integrator.f

    _f = if isinplace(f)
        (resid, u_next, p) -> f(resid, u_next, p.u, p.p, p.t_next)
    else
        (u_next, p) -> f(u_next, p.u, p.p, p.t_next)
    end
    u_len = isnothing(integrator.u) ? 0 : length(integrator.u)
    nlls = !isnothing(f.resid_prototype) && (length(f.resid_prototype) != u_len)

    prob = if nlls
        NonlinearLeastSquaresProblem{isinplace(f)}(
            NonlinearFunction(_f; resid_prototype = f.resid_prototype),
            cache.state.u, cache.state)
    else
        NonlinearProblem{isinplace(f)}(_f, cache.state.u, cache.state)
    end
    cache.prob = prob
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
            (resid, u_next, p) -> f(resid, u_next, p.u, p.p, p.t_next)
        else
            (u_next, p) -> f(u_next, p.u, p.p, p.t_next)
        end

        nlls = !isnothing(f.resid_prototype) &&
               (length(f.resid_prototype) != length(integrator.u))
        prob = if nlls
            NonlinearLeastSquaresProblem{isinplace(f)}(
                NonlinearFunction(_f; resid_prototype = f.resid_prototype), u, initstate)
        else
            NonlinearProblem{isinplace(f)}(_f, u, initstate)
        end
        sol = solve(prob, SimpleNewtonRaphson())
        if sol.retcode == ReturnCode.Success
            integrator.u = sol
        else
            integrator.sol = SciMLBase.solution_new_retcode(integrator.sol, ReturnCode.InitialFailure)
        end
    end
end
