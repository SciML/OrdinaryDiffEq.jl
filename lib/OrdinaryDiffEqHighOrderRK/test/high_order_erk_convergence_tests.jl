# This definitely needs cleaning
using OrdinaryDiffEqHighOrderRK, ODEProblemLibrary, DiffEqDevTools
using Test, Random
Random.seed!(100)

dts5 = 1 .// 2 .^ (3:-1:1)
testTol = 0.2

@testset "Explicit Solver Convergence Tests ($(["out-of-place", "in-place"][i]))" for i in 1:2
    prob = (ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear)[i]

    sim3 = test_convergence(dts5, prob, PFRK87())
    @test sim3.𝒪est[:l∞]≈8.4 atol=0.2
    
    # Test for FW15Stage10 (10th order method)
    dts10 = 1 .// 2 .^ (2:-1:0)  # Fewer steps needed for higher order
    sim4 = test_convergence(dts10, prob, FW15Stage10(), adaptive=false)
    @test sim4.𝒪est[:l∞]≈10.0 atol=1.0  # Order 10 method
end
