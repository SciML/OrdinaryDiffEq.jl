# This definitely needs cleaning
using OrdinaryDiffEqPRK, ODEProblemLibrary, DiffEqDevTools
using Test, Random
Random.seed!(100)

dts2 = 1 .// 2 .^ (7:-1:3)
testTol = 0.2

@testset "Explicit Solver Convergence Tests ($(["out-of-place", "in-place"][i]))" for i in 1:2
    prob = (
        ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear,
    )[i]
    sim3 = test_convergence(dts2, prob, KuttaPRK2p5(threading = true))
    @test sim3.𝒪est[:l∞] ≈ 5 atol = testTol

    sim3 = test_convergence(dts2, prob, KuttaPRK2p5(threading = false))
    @test sim3.𝒪est[:l∞] ≈ 5 atol = testTol
end
