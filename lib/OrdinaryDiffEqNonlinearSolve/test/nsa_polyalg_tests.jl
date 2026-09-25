using OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using NonlinearSolve: RobustMultiNewton, NewtonRaphson
using ADTypes, LinearAlgebra, SciMLBase
using Test

# Regression test for #3859: a NonlinearSolve polyalgorithm (e.g. RobustMultiNewton)
# used as the NonlinearSolveAlg inner solver used to crash in `initialize!` with
# `type NonlinearSolvePolyAlgorithmCache has no field u`, because a polyalg cache keeps
# `u`/`fu` on its active branch rather than at top level; newton.jl now reads them
# through `NonlinearSolveBase.get_u`/`get_fu`. A second bug lived upstream: the polyalg
# cache's `reinit!` did not clear `force_stop`/`retcode`, so once any branch converged,
# every later `step!` short-circuited on `not_terminated`, silently reusing a stale `u`
# and diverging (fixed in NonlinearSolveBase 2.39). This test exercises both: it would
# throw without the accessors and return wrong values without the `reinit!` reset.

const A = [-1000.0 1.0; 1.0 -2.0]
function lin!(du, u, p, t)
    mul!(du, A, u)
    return nothing
end
prob = ODEProblem(lin!, [1.0, 1.0], (0.0, 0.5))
uexact = exp(A * 0.5) * prob.u0

poly = NonlinearSolveAlg(RobustMultiNewton(; autodiff = AutoForwardDiff()))
single = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoForwardDiff()))

# Each run is checked against the analytic solution rather than the polyalg run against
# the NewtonRaphson one: the outer iteration takes exactly one inner `step!` per
# iteration, and RobustMultiNewton's active trust-region branch damps that step, so the
# two inner solvers legitimately produce different step sequences (and step counts).
# The stale-iterate bug produced O(1) errors, far above this bound.
@testset "polyalg inner solver solves as accurately as single-algorithm" begin
    for ODEAlg in (TRBDF2, FBDF)
        ref = solve(prob, ODEAlg(nlsolve = single); reltol = 1.0e-8, abstol = 1.0e-10)
        sol = solve(prob, ODEAlg(nlsolve = poly); reltol = 1.0e-8, abstol = 1.0e-10)
        @test SciMLBase.successful_retcode(ref)
        @test SciMLBase.successful_retcode(sol)
        @test norm(ref.u[end] .- uexact) / norm(uexact) < 1.0e-4
        @test norm(sol.u[end] .- uexact) / norm(uexact) < 1.0e-4
    end
end
