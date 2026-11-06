using OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using NonlinearSolve: DFSane, NewtonRaphson, Broyden
using LineSearch: BackTracking
using ADTypes, LinearAlgebra, SciMLBase
using Test

# The outer convergence test used to infer "close to the root" from the displacement the inner
# solver produced, and never from the residual. That inference is only sound for a step that
# solves the Newton system, where the displacement is `W⁻¹F`. An inner solver whose line search
# collapses leaves the iterate unmoved — either stopping the cache with a failure retcode that
# `compute_step!` only looked for on *entry*, or returning a step orders of magnitude below
# roundoff with the residual bit-for-bit unchanged. Either way `ndz` is tiny for a reason that
# has nothing to do with the distance to the root, `iter == 1 && ndz < 1e-5` fires, and the
# stage is accepted unconverged.
#
# `DFSane` and a back-tracking `Broyden` genuinely cannot solve these stages. The property
# under test is only that the failure reaches the step controller instead of being dressed up
# as a converged solve: on unfixed code ROBER + FBDF + DFSane runs 21 accepted steps and
# returns the initial condition with `ReturnCode.Success`.

function rober!(du, u, p, t)
    y1, y2, y3 = u
    k1, k2, k3 = p
    du[1] = -k1 * y1 + k3 * y2 * y3
    du[2] = k1 * y1 - k2 * y2^2 - k3 * y2 * y3
    du[3] = k2 * y2^2
    return nothing
end
const rober_prob = ODEProblem(
    rober!, [1.0, 0.0, 0.0], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4)
)

function hires!(du, u, p, t)
    y1, y2, y3, y4, y5, y6, y7, y8 = u
    du[1] = -1.71y1 + 0.43y2 + 8.32y3 + 0.0007
    du[2] = 1.71y1 - 8.75y2
    du[3] = -10.03y3 + 0.43y4 + 0.035y5
    du[4] = 8.32y2 + 1.71y3 - 1.12y4
    du[5] = -1.745y5 + 0.43y6 + 0.43y7
    du[6] = -280.0y6 * y8 + 0.69y4 + 1.71y5 - 0.43y6 + 0.69y7
    du[7] = 280.0y6 * y8 - 1.81y7
    du[8] = -280.0y6 * y8 + 1.81y7
    return nothing
end
const hires_prob = ODEProblem(
    hires!, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0057], (0.0, 321.8122)
)

reference(prob) = solve(prob, Rodas5P(); reltol = 1.0e-13, abstol = 1.0e-15).u[end]
relerr(sol, ref) = norm(sol.u[end] .- ref) / norm(ref)

@testset "an unconverged stage is never reported as a successful solve" begin
    for (prob, ODEAlg, inner) in (
            (rober_prob, FBDF, () -> DFSane()),
            (rober_prob, KenCarp4, () -> DFSane()),
            (hires_prob, FBDF, () -> DFSane()),
            (hires_prob, FBDF, () -> Broyden(; linesearch = BackTracking())),
            (hires_prob, KenCarp4, () -> Broyden(; linesearch = BackTracking())),
        )
        ref = reference(prob)
        sol = solve(
            prob, ODEAlg(nlsolve = NonlinearSolveAlg(inner()));
            reltol = 1.0e-6, abstol = 1.0e-9, maxiters = 10^5
        )
        @test !SciMLBase.successful_retcode(sol) || relerr(sol, ref) < 1.0e-3
    end
end

@testset "a solvable stage still converges" begin
    for prob in (rober_prob, hires_prob), ODEAlg in (FBDF, KenCarp4)
        ref = reference(prob)
        sol = solve(
            prob,
            ODEAlg(nlsolve = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoForwardDiff())));
            reltol = 1.0e-6, abstol = 1.0e-9
        )
        @test SciMLBase.successful_retcode(sol)
        @test relerr(sol, ref) < 1.0e-4
    end
end
