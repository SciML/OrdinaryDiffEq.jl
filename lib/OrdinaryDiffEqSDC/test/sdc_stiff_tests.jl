# The convergence gate runs on non-stiff problems, which is the regime where the
# choice of QΔ does not matter. Here it is the only thing that matters: QΔ sets
# the stiff-limit iteration matrix `K_S = I - QΔ⁻¹Q`, so LU (whose `K_S` is
# nilpotent by construction) should stay accurate as λ → -∞ while the Euler
# sweepers degrade and the explicit ones diverge.

using OrdinaryDiffEqSDC
using SciMLBase
using Test

# Prothero-Robinson has the exact solution `cos t` for every λ, so stiffness can
# be dialled without changing the answer.
prothero_robinson(λ) = ODEProblem(
    (u, p, t) -> λ * (u - cos(t)) - sin(t), 1.0, (0.0, 1.0)
)

function pr_error(sweeper, λ, K; M = 3, nsteps = 10)
    alg = SDC(num_nodes = M, num_sweeps = K, sweeper = sweeper)
    sol = solve(
        prothero_robinson(λ), alg; dt = 1 / nsteps, adaptive = false,
        save_everystep = true, dense = false
    )
    return maximum(abs(u - cos(t)) for (u, t) in zip(sol.u, sol.t))
end

@testset "LU sweeper stays accurate under stiffness" begin
    # The LU trick makes the stiff-limit iteration matrix nilpotent, so the
    # error should be near enough independent of λ.
    for λ in (-1.0e2, -1.0e4, -1.0e6)
        @test pr_error(SDCSweeper.LU, λ, 3) < 1.0e-2
        @test pr_error(SDCSweeper.LU, λ, 6) < 1.0e-5
    end
    # Same step size, four more orders of magnitude of stiffness.
    @test pr_error(SDCSweeper.LU, -1.0e6, 3) < 10 * pr_error(SDCSweeper.LU, -1.0e2, 3)
end

@testset "sweeper ranking under stiffness" begin
    λ = -1.0e4
    be = pr_error(SDCSweeper.BE, λ, 3)
    lu = pr_error(SDCSweeper.LU, λ, 3)
    bepar = pr_error(SDCSweeper.BEpar, λ, 3)

    # This gap is the whole point of the LU trick and is invisible on a
    # non-stiff problem.
    @test lu < be / 100

    # The naive diagonal sweeper is worse than the lower-triangular one, which is
    # the cost parallel-across-the-nodes SDC has to buy back with better
    # coefficients rather than with the parallelism alone.
    @test bepar > be

    # Explicit sweeps cannot integrate a stiff problem at this step size.
    @test !(pr_error(SDCSweeper.FE, λ, 3) < 1)
    @test !(pr_error(SDCSweeper.Trapezoid, λ, 3) < 1)
end

@testset "stiff decay is resolved, not amplified" begin
    # dt·|λ| = 10⁶: nothing here is resolved, so an implicit method must damp it
    # and an explicit one must blow up.
    decay = ODEProblem((u, p, t) -> -1.0e6 * u, 1.0, (0.0, 10.0))
    for sweeper in (SDCSweeper.BE, SDCSweeper.LU)
        sol = solve(
            decay, SDC(num_nodes = 3, num_sweeps = 3, sweeper = sweeper);
            dt = 1.0, adaptive = false
        )
        @test all(abs.(sol.u[2:end]) .< abs.(sol.u[1:(end - 1)]))
        @test abs(sol.u[end]) < 1.0e-2
    end
    sol = solve(
        decay, SDC(num_nodes = 3, num_sweeps = 3, sweeper = SDCSweeper.FE);
        dt = 1.0, adaptive = false
    )
    @test !isfinite(sol.u[end]) || abs(sol.u[end]) > 1
end

@testset "stiff system, in-place, autodiff Jacobian" begin
    function stiff!(du, u, p, t)
        du[1] = -1.0e4 * u[1]
        du[2] = -u[2]
        return nothing
    end
    prob = ODEProblem(stiff!, [1.0, 1.0], (0.0, 1.0))
    sol = solve(
        prob, SDC(num_nodes = 3, num_sweeps = 4, sweeper = SDCSweeper.LU);
        dt = 0.1, adaptive = false, save_everystep = true, dense = false
    )
    # The fast mode is unresolved, so it only has to be damped. The slow mode
    # must still reach full order, which is what fails when the stiff component
    # contaminates the sweep.
    fast = [abs(u[1]) for u in sol.u]
    @test all(fast[2:end] .< fast[1:(end - 1)])
    @test fast[end] < 1.0e-8
    @test maximum(abs(u[2] - exp(-t)) for (u, t) in zip(sol.u, sol.t)) < 1.0e-6

    # LU damps the unresolved mode far harder per step than backward Euler.
    be = solve(
        prob, SDC(num_nodes = 3, num_sweeps = 4, sweeper = SDCSweeper.BE);
        dt = 0.1, adaptive = false, save_everystep = true, dense = false
    )
    @test fast[2] < abs(be.u[2][1]) / 10
end
