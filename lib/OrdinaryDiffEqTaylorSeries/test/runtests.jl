using OrdinaryDiffEqTaylorSeries, ODEProblemLibrary, DiffEqDevTools
using Test

@testset "Taylor2 Convergence Tests" begin
    # Test convergence
    dts = 2.0 .^ (-8:-4)
    testTol = 0.2
    sim = test_convergence(dts, prob_ode_linear, ExplicitTaylor2())
    @test sim.𝒪est[:final]≈2 atol=testTol
    sim = test_convergence(dts, prob_ode_2Dlinear, ExplicitTaylor2())
    @test sim.𝒪est[:final]≈2 atol=testTol
end

@testset "Taylor Convergence Tests" begin
    # Test convergence
    dts = 2.0 .^ (-8:-4)
    testTol = 0.2
    for N in 3:4
        alg = ExplicitTaylor(order = Val(N))
        sim = test_convergence(dts, prob_ode_linear, alg)
        @test sim.𝒪est[:final]≈N atol=testTol
        sim = test_convergence(dts, prob_ode_2Dlinear, alg)
        @test sim.𝒪est[:final]≈N atol=testTol
    end
end

@testset "Taylor Adaptive time-step Tests" begin
    sol = solve(prob_ode_linear, ExplicitTaylor(order = Val(4)))
    @test length(sol) < 20
end

@testset "Taylor Adaptive time-step Adaptive order Tests" begin
    sol = solve(prob_ode_linear, ExplicitTaylorAdaptiveOrder())
    @test length(sol) < 20
end
