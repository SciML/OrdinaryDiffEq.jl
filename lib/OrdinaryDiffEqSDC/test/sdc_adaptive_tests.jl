using OrdinaryDiffEqSDC
using SciMLBase
using Test

rotation!(du, u, p, t) = (du[1] = -u[2]; du[2] = u[1]; nothing)
const ROTATION = ODEProblem(rotation!, [1.0, 0.0], (0.0, 2 * π))
rotation_exact(t) = [cos(t), sin(t)]

final_error(sol) = maximum(abs.(sol.u[end] .- rotation_exact(sol.t[end])))

@testset "SDC adaptive order bookkeeping" begin
    # The embedded solution is the step update one sweep earlier.
    for K in 1:5
        alg = SDC(num_nodes = 4, num_sweeps = K)
        @test OrdinaryDiffEqSDC.alg_adaptive_order(alg) ==
            max(1, OrdinaryDiffEqSDC.alg_order(alg) - 1)
    end
end

@testset "SDC honours the requested tolerance" begin
    alg = SDC(num_nodes = 3, num_sweeps = 3)
    errors = Float64[]
    steps = Int[]
    for tol in (1.0e-4, 1.0e-6, 1.0e-8)
        sol = solve(ROTATION, alg; abstol = tol, reltol = tol)
        @test SciMLBase.successful_retcode(sol)
        push!(errors, final_error(sol))
        push!(steps, length(sol.t) - 1)
        # A factor of 100 leaves room for the estimator being a local one while
        # the measurement is a global error.
        @test errors[end] < 100 * tol
    end
    @test issorted(errors, rev = true)
    @test issorted(steps)
end

@testset "SDC adaptive reuses Jacobians across sweeps" begin
    # Fixed step size makes the library refactorise on every nonlinear solve
    # (`do_newJW`: `!integrator.opts.adaptive && return true, true`), which for
    # SDC is one per node per sweep; adaptive steps reuse instead.
    alg = SDC(num_nodes = 3, num_sweeps = 3)
    adaptive = solve(ROTATION, alg; abstol = 1.0e-8, reltol = 1.0e-8)
    nsteps = length(adaptive.t) - 1
    fixed = solve(ROTATION, alg; dt = (2 * π) / nsteps, adaptive = false)

    @test fixed.stats.njacs == nsteps * alg.num_nodes * alg.num_sweeps
    @test adaptive.stats.njacs < nsteps
    @test adaptive.stats.nw < fixed.stats.nw / 10
end

@testset "SDC adaptive out-of-place" begin
    prob = ODEProblem((u, p, t) -> [-u[2], u[1]], [1.0, 0.0], (0.0, 2 * π))
    sol = solve(prob, SDC(num_nodes = 3, num_sweeps = 3); abstol = 1.0e-8, reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
    @test maximum(abs.(sol.u[end] .- rotation_exact(sol.t[end]))) < 1.0e-6
end

@testset "SDC adaptive on a stiff problem needs a converged iteration" begin
    # The estimate is the difference between two sweeps, so it measures how far
    # the iteration is from the collocation solution rather than how far that
    # solution is from the truth. On a stiff problem the iteration converges
    # slowly, so too few sweeps make the estimate — not the accuracy — set the
    # step size. Enough sweeps and the step size is accuracy limited again.
    prothero_robinson = ODEProblem(
        (u, p, t) -> -1.0e4 * (u - cos(t)) - sin(t), 1.0, (0.0, 1.0)
    )
    solve_with(K) = solve(
        prothero_robinson,
        SDC(num_nodes = 3, num_sweeps = K, sweeper = SDCSweeper.LU);
        abstol = 1.0e-8, reltol = 1.0e-8
    )
    few, many = solve_with(3), solve_with(8)
    for sol in (few, many)
        @test SciMLBase.successful_retcode(sol)
        @test maximum(abs(u - cos(t)) for (u, t) in zip(sol.u, sol.t)) < 1.0e-6
    end
    @test length(many.t) < length(few.t) / 10
    @test many.stats.nsolve < few.stats.nsolve / 10
end
