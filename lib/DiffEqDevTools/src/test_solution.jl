"""
TestSolution

"""
mutable struct TestSolution{T, N, hasinterp, tType, uType, iType} <:
               AbstractTimeseriesSolution{T, N, uType}
    t::tType
    u::uType
    interp::iType
    dense::Bool
    retcode::ReturnCode.T
end
(T::TestSolution)(t) = T.interp(t)
function TestSolution(t, u)
    T = eltype(eltype(u))
    N = length((size(u[1])..., length(u)))
    TestSolution{T, N, false, typeof(t), typeof(u), Nothing}(t, u, nothing, false,
        ReturnCode.Success)
end
function TestSolution(t, u, interp)
    T = eltype(eltype(u))
    N = length((size(u[1])..., length(u)))
    TestSolution{T, N, true, typeof(t), typeof(u), typeof(interp)}(t, u, interp, true,
        ReturnCode.Success)
end
function TestSolution(interp::DESolution)
    TestSolution{Nothing, 0, true, Nothing, Nothing, typeof(interp)}(nothing, nothing,
        interp, true,
        ReturnCode.Success)
end
function hasinterp(::TestSolution{
        T, N, hi, tType, uType, iType}) where {T, N, hi, tType,
        uType, iType}
    hi
end
"""
`appxtrue(sol::AbstractODESolution,sol2::TestSolution)`

Uses the interpolant from the higher order solution sol2 to approximate
errors for sol. If sol2 has no interpolant, only the final error is
calculated.
"""
function appxtrue(sol::AbstractODESolution, sol2::TestSolution)
    if sol2.u == nothing && hasinterp(sol2)
        _sol = TestSolution(sol.t, sol2(sol.t).u, sol2)
    else
        _sol = sol2
    end

    errors = Dict(:final => recursive_mean(abs.(sol.u[end] - _sol.u[end])))
    if _sol.dense
        timeseries_analytic = _sol(sol.t)
        errors[:l∞] = maximum(vecvecapply((x) -> abs.(x), sol - timeseries_analytic))
        errors[:l2] = sqrt(recursive_mean(vecvecapply((x) -> float(x) .^ 2,
            sol - timeseries_analytic)))
        densetimes = collect(range(sol.t[1], stop = sol.t[end], length = 100))
        interp_u = sol(densetimes)
        interp_analytic = _sol(densetimes)
        interp_errors = Dict(
            :L∞ => maximum(vecvecapply((x) -> abs.(x),
                interp_u - interp_analytic)),
            :L2 => sqrt(recursive_mean(vecvecapply((x) -> float(x) .^ 2,
                interp_u -
                interp_analytic))))
        errors = merge(errors, interp_errors)
    else
        timeseries_analytic = sol2.u
        if sol.t == sol2.t
            errors[:l∞] = maximum(vecvecapply((x) -> abs.(x), sol - timeseries_analytic))
            errors[:l2] = sqrt(recursive_mean(vecvecapply((x) -> float(x) .^ 2,
                sol - timeseries_analytic)))
        end
    end
    DiffEqBase.build_solution(sol, timeseries_analytic, errors)
end

"""
`appxtrue(sol::AbstractODESolution,sol2::AbstractODESolution)`

Uses the interpolant from the higher order solution sol2 to approximate
errors for sol. If sol2 has no interpolant, only the final error is
calculated.
"""
function appxtrue(sol::AbstractODESolution, sol2::AbstractODESolution;
        timeseries_errors = sol2.dense, dense_errors = sol2.dense)
    errors = Dict(:final => recursive_mean(abs.(sol.u[end] - sol2.u[end])))
    if sol2.dense
        timeseries_analytic = sol2(sol.t)
        errors[:l∞] = maximum(vecvecapply((x) -> abs.(x), sol - timeseries_analytic))
        errors[:l2] = sqrt(recursive_mean(vecvecapply((x) -> float(x) .^ 2,
            sol - timeseries_analytic)))
        if dense_errors
            densetimes = collect(range(sol.t[1], stop = sol.t[end], length = 100))
            interp_u = sol(densetimes)
            interp_analytic = sol2(densetimes)
            interp_errors = Dict(
                :L∞ => maximum(vecvecapply((x) -> abs.(x),
                    interp_u - interp_analytic)),
                :L2 => sqrt(recursive_mean(vecvecapply((x) -> float(x) .^ 2,
                    interp_u -
                    interp_analytic))))
            errors = merge(errors, interp_errors)
        end
    else
        timeseries_analytic = sol2.u
        if timeseries_errors && sol.t == sol2.t
            errors[:l∞] = maximum(vecvecapply((x) -> abs.(x), sol - timeseries_analytic))
            errors[:l2] = sqrt(recursive_mean(vecvecapply((x) -> float(x) .^ 2,
                sol - timeseries_analytic)))
        end
    end
    DiffEqBase.build_solution(sol, timeseries_analytic, errors)
end

function appxtrue(sim::EnsembleSolution, appx_setup; kwargs...)
    _new_sols = Vector{DESolution}(length(sim.u))
    for i in eachindex(sim)
        prob = sim[i].prob
        prob2 = SDEProblem(prob.f, prob.g, prob.u0, prob.tspan,
            noise = NoiseWrapper(sim.u[i].W))
        true_sol = solve(prob2, appx_setup[:alg]; appx_setup...)
        _new_sols[i] = appxtrue(sim.u[i], true_sol)
    end
    new_sols = convert(Vector{typeof(_new_sols[1])}, _new_sols)
    calculate_ensemble_errors(new_sols; converged = sim.converged,
        elapsedTime = sim.elapsedTime, kwargs...)
end
