# This definitely needs cleaning
using OrdinaryDiffEqPDIRK, ODEProblemLibrary, DiffEqDevTools
using Test, Random
Random.seed!(100)

## Convergence Testing
testTol = 0.2

@testset "Implicit Solver Convergence Tests ($(["out-of-place", "in-place"][i]))" for i in 1:2
    dts = 1 .// 2 .^ (9:-1:5)
    prob = (
        ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear,
    )[i]
    sim18 = test_convergence(dts, prob, PDIRK44())
    @test sim18.𝒪est[:final] ≈ 4 atol = testTol

    sim182 = test_convergence(dts, prob, PDIRK44(; threading = false))
    @test sim182.𝒪est[:final] ≈ 4 atol = testTol
end
