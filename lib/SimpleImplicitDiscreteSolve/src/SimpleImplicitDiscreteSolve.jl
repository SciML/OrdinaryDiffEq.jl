module SimpleImplicitDiscreteSolve

using SciMLBase
using SimpleNonlinearSolve
using Reexport
using StaticArrays
@reexport using DiffEqBase

"""
    SimpleIDSolve()

Simple solver for `ImplicitDiscreteSystems`. Uses `SimpleNewtonRaphson` to solve for the next state at every timestep.
"""
struct SimpleIDSolve <: SciMLBase.AbstractODEAlgorithm end

function DiffEqBase.__init(prob::ImplicitDiscreteProblem, alg::SimpleIDSolve; dt = 1)
    u0 = prob.u0
    p = prob.p
    f = prob.f
    t = prob.tspan[1]

    nlf = isinplace(f) ? (out, u, p) -> f(out, u, u0, p, t) : (u, p) -> f(u, u0, p, t)
    prob = NonlinearProblem{isinplace(f)}(nlf, u0, p)
    sol = solve(prob, SimpleNewtonRaphson())
    sol, (sol.retcode != ReturnCode.Success)
end

function DiffEqBase.solve(prob::ImplicitDiscreteProblem, alg::SimpleIDSolve;
        dt = 1,
        save_everystep = true,
        save_start = true,
        adaptive = false,
        dense = false,
        save_end = true,
        kwargs...)
    @assert !adaptive
    @assert !dense
    (initsol, initfail) = DiffEqBase.__init(prob, alg; dt)
    if initfail
        sol = DiffEqBase.build_solution(prob, alg, prob.tspan[1], u0, k = nothing,
            stats = nothing, calculate_error = false)
        return SciMLBase.solution_new_retcode(sol, ReturnCode.InitialFailure)
    end

    u0 = initsol.u
    tspan = prob.tspan
    f = prob.f
    p = prob.p
    t = tspan[1]
    tf = prob.tspan[2]
    ts = tspan[1]:dt:tspan[2]

    l = save_everystep ? length(ts) - 1 : 1
    save_start && (l = l + 1)
    u0type = typeof(u0)
    us = u0type <: StaticArray ? MVector{l, u0type}(undef) : Vector{u0type}(undef, l)

    if save_start
        us[1] = u0
    end

    u = u0
    convfail = false
    for i in 2:length(ts)
        uprev = u
        t = ts[i]
        nlf = isinplace(f) ? (out, u, p) -> f(out, u, uprev, p, t) :
              (u, p) -> f(u, uprev, p, t)
        nlprob = NonlinearProblem{isinplace(f)}(nlf, uprev, p)
        nlsol = solve(nlprob, SimpleNewtonRaphson())
        u = nlsol.u
        save_everystep && (us[i] = u)
        convfail = (nlsol.retcode != ReturnCode.Success)

        if convfail
            sol = DiffEqBase.build_solution(prob, alg, ts[1:i], us[1:i], k = nothing,
                stats = nothing, calculate_error = false)
            sol = SciMLBase.solution_new_retcode(sol, ReturnCode.ConvergenceFailure)
            return sol
        end
    end

    !save_everystep && save_end && (us[end] = u)
    sol = DiffEqBase.build_solution(prob, alg, ts, us,
        k = nothing, stats = nothing,
        calculate_error = false)

    DiffEqBase.has_analytic(prob.f) &&
        DiffEqBase.calculate_solution_errors!(
            sol; timeseries_errors = true, dense_errors = false)
    sol
end

export SimpleIDSolve

end
