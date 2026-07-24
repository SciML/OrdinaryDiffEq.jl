using GlobalDiffEq, OrdinaryDiffEqTsit5, LinearAlgebra
using RecursiveArrayTools: ArrayPartition
using Test
import SciMLBase

# Unstable Prince42 problem: y' = y - sin(t) + cos(t), y(0) = 0, exact y = sin(t).
# Perturbations grow like e^t, making it the classic global-error-estimation test.
f_prince!(du, u, p, t) = (du[1] = u[1] - sin(t) + cos(t); nothing)
f_prince(u, p, t) = [u[1] - sin(t) + cos(t)]
const prince_tspan = (0.0, 2.0)
prince_exact(t) = sin(t)

function glee_convergence(prob, alg, dts)
    errs = Float64[]
    est_errs = Float64[]
    for dt in dts
        sol = solve(prob, alg; dt = dt, adaptive = false)
        err = prince_exact(prince_tspan[2]) - sol.u[end].x[1][1]
        est = global_error_estimate(sol, length(sol.u))[1]
        push!(errs, abs(err))
        push!(est_errs, abs(est - err))
    end
    sol_order = log2(errs[end - 1] / errs[end])
    est_order = log2(est_errs[end - 1] / est_errs[end])
    return sol_order, est_order
end

@testset "GLEE convergence and estimate accuracy" begin
    dts = 2.0 .^ -(5:8)
    for iip in (true, false)
        prob = if iip
            ODEProblem(f_prince!, [0.0], prince_tspan)
        else
            ODEProblem(f_prince, [0.0], prince_tspan)
        end
        for (alg, p) in (
                (GLEE23(), 2), (GLEE24(), 2), (GLEE35(), 3),
            )
            sol_order, est_order = glee_convergence(prob, alg, dts)
            @test sol_order ≈ p atol = 0.15
            # ε is an asymptotically correct global error estimate: it converges
            # to the true error one order faster than the solution converges
            @test est_order ≈ p + 1 atol = 0.2
        end
    end
end

@testset "GLEE estimate tracks the true error" begin
    prob = ODEProblem(f_prince!, [0.0], prince_tspan)
    for alg in (GLEE23(), GLEE24(), GLEE35())
        sol = solve(prob, alg; dt = 1 / 128, adaptive = false)
        err = prince_exact(prince_tspan[2]) - sol.u[end].x[1][1]
        est = global_error_estimate(sol)[end][1]
        @test est / err ≈ 1 rtol = 0.1
    end
end

@testset "GLEE adaptive stepping" begin
    prob = ODEProblem(f_prince!, [0.0], prince_tspan)
    for alg in (GLEE24(), GLEE35())
        sol = solve(prob, alg; abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        err = prince_exact(prince_tspan[2]) - sol.u[end].x[1][1]
        est = global_error_estimate(sol)[end][1]
        @test isfinite(est)
        @test est / err ≈ 1 rtol = 0.6
    end
end

@testset "GLEE integrator interface" begin
    prob = ODEProblem(f_prince!, [0.0], prince_tspan)
    integrator = init(prob, GLEE24(); dt = 1 / 64, adaptive = false)
    @test integrator.u isa ArrayPartition
    solve!(integrator)
    sol = integrator.sol
    @test SciMLBase.successful_retcode(sol)
    err = prince_exact(prince_tspan[2]) - sol.u[end].x[1][1]
    @test global_error_estimate(sol)[end][1] / err ≈ 1 rtol = 0.1
end

@testset "GLEE argument validation" begin
    prob = ODEProblem(f_prince!, [0.0], prince_tspan)
    partitioned_prob = ODEProblem(
        (du, u, p, t) -> nothing, ArrayPartition([0.0], [0.0]), prince_tspan
    )
    @test_throws ArgumentError solve(partitioned_prob, GLEE24(); dt = 0.1)

    mm_prob = ODEProblem(
        ODEFunction(f_prince!; mass_matrix = [2.0;;]), [0.0], prince_tspan
    )
    @test_throws ArgumentError solve(mm_prob, GLEE24(); dt = 0.1)

    sol = solve(prob, Tsit5())
    @test_throws ArgumentError global_error_estimate(sol, 1)

    @test SciMLBase.alg_order(GLEE23()) == 2
    @test SciMLBase.alg_order(GLEE24()) == 2
    @test SciMLBase.alg_order(GLEE35()) == 3
end
