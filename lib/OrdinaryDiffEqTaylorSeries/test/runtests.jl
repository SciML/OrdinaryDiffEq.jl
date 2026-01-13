using OrdinaryDiffEqTaylorSeries, ODEProblemLibrary, DiffEqDevTools
using Test
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @testset "Taylor2 Convergence Tests" begin
        # Test convergence
        dts = 2.0 .^ (-8:-4)
        testTol = 0.2
        sim = test_convergence(dts, prob_ode_linear, ExplicitTaylor2())
        @test sim.𝒪est[:final] ≈ 2 atol = testTol
        sim = test_convergence(dts, prob_ode_2Dlinear, ExplicitTaylor2())
        @test sim.𝒪est[:final] ≈ 2 atol = testTol
    end

    @testset "TaylorN Convergence Tests" begin
        # Test convergence
        dts = 2.0 .^ (-8:-4)
        testTol = 0.2
        for N in 3:4
            alg = ExplicitTaylor(order = Val(N))
            sim = test_convergence(dts, prob_ode_linear, alg)
            @test sim.𝒪est[:final] ≈ N atol = testTol
            sim = test_convergence(dts, prob_ode_2Dlinear, alg)
            @test sim.𝒪est[:final] ≈ N atol = testTol
        end
    end

    @testset "TaylorN Adaptive Tests" begin
        sol = solve(prob_ode_linear, ExplicitTaylor(order = Val(2)))
        @test length(sol) < 20
        @test SciMLBase.successful_retcode(sol)
    end
end

# Run QA tests (JET, Aqua)
if TEST_GROUP != "FUNCTIONAL"
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
