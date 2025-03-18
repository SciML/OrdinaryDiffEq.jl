# Remake the nonlinear problem, then update
function perform_step!(integrator, cache::SimpleIDSolveCache, repeat_step = false)
    @unpack alg, u, uprev, dt, t, f, p = integrator
    @unpack state, prob = cache
    state.u .= uprev
    state.t_next = t
    prob = remake(prob, p = state)

    u = solve(prob, SimpleNewtonRaphson())
    integrator.sol = SciMLBase.solution_new_retcode(integrator.sol, u.retcode)
    integrator.u = u
end

function initialize!(integrator, cache::SimpleIDSolveCache)
    integrator.u isa AbstractVector && (cache.state.u .= integrator.u)
    cache.state.p = integrator.p
    cache.state.t_next = integrator.t
    f = integrator.f

    _f = if isinplace(f)
        (resid, u_next, p) -> f(resid, u_next, p.u, p.p, p.t_next)
    else
        (u_next, p) -> f(u_next, p.u, p.p, p.t_next)
    end

    prob = if isinplace(f)
        NonlinearProblem{true}(_f, cache.state.u, cache.state)
    else
        NonlinearProblem{false}(_f, cache.state.u, cache.state)
    end
    cache.prob = prob
end

function _initialize_dae!(integrator, prob::ImplicitDiscreteProblem,
        alg::DefaultInit, x::Union{Val{true}, Val{false}})
    isnothing(prob.u0) && return
    atol = one(eltype(u0)) * 1e-12
    if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(integrator, prob,
                         OverrideInit(atol), x)
    else
        @unpack u, p, t, f = integrator
        initstate = ImplicitDiscreteState(u, p, t)

        _f = if isinplace(f)
            (resid, u_next, p) -> f(resid, u_next, p.u, p.p, p.t_next)
        else
            (u_next, p) -> f(u_next, p.u, p.p, p.t_next)
        end
        prob = NonlinearProblem{isinplace(f)}(_f, u, initstate)
        sol = solve(prob, SimpleNewtonRaphson())
        integrator.u = sol
    end
end
