# This definitely needs cleaning
using OrdinaryDiffEq, ODEProblemLibrary, DiffEqDevTools
using Test, Random
Random.seed!(100)

testTol = 0.2

@testset "Implicit Solver Convergence Tests ($(["out-of-place", "in-place"][i]))" for i in 1:2
    prob = (ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear)[i]

    sim = test_convergence(dts, prob, ABDF2())
    @test sim.𝒪est[:final]≈2 atol=testTol
    @test sim.𝒪est[:l2]≈2 atol=testTol
    @test sim.𝒪est[:l∞]≈2 atol=testTol

    sim = test_convergence(dts, prob, ABDF2(nlsolve = NLFunctional()))
    @test sim.𝒪est[:final]≈2 atol=testTol
    @test sim.𝒪est[:l2]≈2 atol=testTol
    @test sim.𝒪est[:l∞]≈2 atol=testTol

    # QBDF
    sim = test_convergence(dts, prob, QBDF1())
    @test sim.𝒪est[:final]≈1 atol=testTol
    @test sim.𝒪est[:l2]≈1 atol=testTol
    @test sim.𝒪est[:l∞]≈1 atol=testTol

    sim = test_convergence(dts, prob, QBDF2())
    @test sim.𝒪est[:final]≈2 atol=testTol
    @test sim.𝒪est[:l2]≈2 atol=testTol
    @test sim.𝒪est[:l∞]≈2 atol=testTol

    # QNDF
    sim = test_convergence(dts, prob, QNDF1())
    @test sim.𝒪est[:final]≈1 atol=testTol
    @test sim.𝒪est[:l2]≈1 atol=testTol
    @test sim.𝒪est[:l∞]≈1 atol=testTol

    sim = test_convergence(dts3, prob, QNDF2())
    @test sim.𝒪est[:final]≈2 atol=testTol
    @test sim.𝒪est[:l2]≈2 atol=testTol
    @test sim.𝒪est[:l∞]≈2 atol=testTol

    sim = test_convergence(dts, prob, QNDF2(nlsolve = NLFunctional()))
    @test sim.𝒪est[:final]≈2 atol=testTol
    @test sim.𝒪est[:l2]≈2 atol=testTol
    @test sim.𝒪est[:l∞]≈2 atol=testTol
    @test_nowarn solve(prob, QNDF())

    # MEBDF2
    sim21 = test_convergence(dts, prob, MEBDF2(extrapolant = :linear))
    @test sim21.𝒪est[:final]≈2 atol=testTol

    sim22 = test_convergence(dts, prob, MEBDF2(nlsolve = NLFunctional()), reltol = 1e-2)
    @test sim22.𝒪est[:final]≈2 atol=testTol

    sim23 = test_convergence(dts, prob, MEBDF2(nlsolve = NLAnderson()), reltol = 1e-2)
    @test sim23.𝒪est[:final]≈2 atol=testTol

    sim24 = test_convergence(
        dts, prob, MEBDF2(nlsolve = NonlinearSolveAlg()), reltol = 1e-2)
    @test sim24.𝒪est[:final]≈2 atol=testTol

    #FBDF
    @test_nowarn solve(prob, FBDF())
end
