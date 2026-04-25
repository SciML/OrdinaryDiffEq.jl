using OrdinaryDiffEqFIRK, DiffEqDevTools, Test, LinearAlgebra
import OrdinaryDiffEqCore
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear

testTol = 0.6

# test orders 2 and 3 for convergence with fixed dt, anything higher is too sensitive to floating point precision
@testset "GaussLegendre: fixed-dt empirical order (s = 2, 3)" begin
    dts = Float64.(1 ./ 2 .^ (5:-1:2))  # 1/32 … 1/4
    for s in 2:3
        alg = GaussLegendre(num_stages = s; maxiters = 100)
        sim = test_convergence(
            dts, prob_ode_linear, alg;
            dense_errors = false,
            abstol = 1.0e-12, reltol = 1.0e-12
        )
        @test sim.𝒪est[:final] ≈ 2 * s atol = testTol
    end
end

# test accuracy of high-order fixed-dt method (s = 4, order 8) empirically
@testset "GaussLegendre: fixed-dt accuracy (s = 4, order 8)" begin
    s = 4
    alg = GaussLegendre(num_stages = s; maxiters = 100)
    sol = solve(
        prob_ode_linear, alg; adaptive = false, dt = 1 // 256,
        abstol = 1.0e-14, reltol = 1.0e-14
    )
    @test sol.retcode == ReturnCode.Success
    exact = prob_ode_linear.u0 * exp(1.01 * (sol.t[end] - sol.t[1]))
    @test isapprox(sol.u[end], exact; rtol = 1.0e-9, atol = 1.0e-12)
end

# test adaptive stepping with tolerance matching
@testset "GaussLegendre: adaptive run matches tolerance" begin
    for s in 2:4
        reltol = 1.0e-6
        abstol = 1.0e-9
        sol = solve(
            prob_ode_linear, GaussLegendre(num_stages = s);
            reltol = reltol, abstol = abstol
        )
        @test sol.retcode == ReturnCode.Success
        exact = prob_ode_linear.u0 * exp(1.01 * (sol.t[end] - sol.t[1]))
        @test isapprox(sol.u[end], exact; rtol = 1.0e-3, atol = 1.0e-6)
    end
end

@testset "GaussLegendre: adaptive controller defaults to PI" begin
    alg = GaussLegendre(num_stages = 3)
    integrator = init(
        prob_ode_linear, alg;
        reltol = 1.0e-6, abstol = 1.0e-9
    )

    controller = integrator.controller_cache.controller
    order = 2 * alg.num_stages
    @test alg.controller === :PI
    @test integrator.controller_cache isa OrdinaryDiffEqCore.PIControllerCache
    @test controller.beta1 ≈ 7 / (10 * order)
    @test controller.beta2 ≈ 2 / (5 * order)
end

# test that Richardson step-doubling tightens step count with tolerance
@testset "GaussLegendre: PI/Richardson tightens step count when tol tightens" begin
    s = 3
    sol_loose = solve(
        prob_ode_linear, GaussLegendre(num_stages = s);
        reltol = 1.0e-3, abstol = 1.0e-6
    )
    sol_tight = solve(
        prob_ode_linear, GaussLegendre(num_stages = s);
        reltol = 1.0e-8, abstol = 1.0e-10
    )
    @test length(sol_tight.t) >= length(sol_loose.t)
end

# test that num_stages = 1 with adaptive throws
@testset "GaussLegendre: num_stages = 1 with adaptive throws" begin
    @test_throws ArgumentError solve(
        prob_ode_linear, GaussLegendre(num_stages = 1);
        reltol = 1.0e-4, abstol = 1.0e-7
    )
end

# test that num_stages = 1 with adaptive = false runs
@testset "GaussLegendre: num_stages = 1 with adaptive = false runs" begin
    sol = solve(
        prob_ode_linear, GaussLegendre(num_stages = 1);
        adaptive = false, dt = 1 // 32
    )
    @test sol.retcode == ReturnCode.Success
end
