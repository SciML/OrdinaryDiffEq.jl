using OrdinaryDiffEqPDIRK, OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NLNewton, NLFunctional, NLAnderson,
    NonlinearSolveAlg, HomotopyNonlinearSolveAlg
using NonlinearSolve: NewtonRaphson
using LinearSolve: KrylovJL_GMRES

stats_counters(stats) = NamedTuple{fieldnames(typeof(stats))}(
    Tuple(getfield(stats, field) for field in fieldnames(typeof(stats)))
)

function stats_rhs!(du, u, p, t)
    @. du = -u * u
    return nothing
end
stats_rhs(u, p, t) = -u .* u

function yielding_stats_rhs!(du, u, p, t)
    yield()
    return stats_rhs!(du, u, p, t)
end
function yielding_stats_rhs(u, p, t)
    yield()
    return stats_rhs(u, p, t)
end

const stats_nonlinear_solvers = (
    NLNewton(), NLFunctional(), NLAnderson(), NonlinearSolveAlg(NewtonRaphson()),
    HomotopyNonlinearSolveAlg(),
)

function threaded_stats_reference()
    references = Dict()
    for rhs in (stats_rhs!, stats_rhs), nl in stats_nonlinear_solvers
        sol = solve(
            ODEProblem(rhs, ones(8), (0.0, 1.0)), PDIRK44(threading = true, nlsolve = nl);
            dt = 0.01, adaptive = false, save_everystep = false
        )
        references[(nameof(rhs), nameof(typeof(nl)))] = stats_counters(sol.stats)
    end
    sol = solve(
        ODEProblem(stats_rhs!, ones(8), (0.0, 1.0)),
        PDIRK44(threading = true, linsolve = KrylovJL_GMRES(), concrete_jac = false);
        dt = 0.01, adaptive = false, save_everystep = false
    )
    references[:krylov] = stats_counters(sol.stats)
    for rhs in (yielding_stats_rhs!, yielding_stats_rhs)
        sol = solve(
            ODEProblem(rhs, ones(8), (0.0, 1.0)),
            PDIRK44(threading = true, nlsolve = NLNewton(always_new = true));
            dt = 0.01, adaptive = false
        )
        references[nameof(rhs)] = stats_counters(sol.stats)
    end
    return references
end
