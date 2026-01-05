# This definitely needs cleaning
using OrdinaryDiffEqBDF, ODEProblemLibrary, DiffEqDevTools, ADTypes, LinearSolve
using OrdinaryDiffEqNonlinearSolve: NLFunctional, NLAnderson, NonlinearSolveAlg
using Test, Random
Random.seed!(100)

testTol = 0.2
dts = 1 .// 2 .^ (9:-1:5)
dts3 = 1 .// 2 .^ (12:-1:7)

if isempty(VERSION.prerelease)
    using Enzyme
end

@testset "Implicit Solver Convergence Tests ($(["out-of-place", "in-place"][i]))" for i in
    1:2

    prob = (
        ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear,
    )[i]

    sim = test_convergence(dts, prob, ABDF2())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:l2] ≈ 2 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 2 atol = testTol

    sim = test_convergence(dts, prob, ABDF2(nlsolve = NLFunctional()))
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:l2] ≈ 2 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 2 atol = testTol

    # QBDF
    sim = test_convergence(dts, prob, QBDF1())
    @test sim.𝒪est[:final] ≈ 1 atol = testTol
    @test sim.𝒪est[:l2] ≈ 1 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 1 atol = testTol

    sim = test_convergence(dts, prob, QBDF2())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:l2] ≈ 2 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 2 atol = testTol

    # QNDF
    sim = test_convergence(dts, prob, QNDF1())
    @test sim.𝒪est[:final] ≈ 1 atol = testTol
    @test sim.𝒪est[:l2] ≈ 1 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 1 atol = testTol

    sim = test_convergence(dts, prob, QNDF1(autodiff = AutoFiniteDiff()))
    @test sim.𝒪est[:final] ≈ 1 atol = testTol
    @test sim.𝒪est[:l2] ≈ 1 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 1 atol = testTol

    if isempty(VERSION.prerelease)
        sim = test_convergence(
            dts,
            prob,
            QNDF1(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward),
                    function_annotation = Enzyme.Const
                )
            )
        )
        @test sim.𝒪est[:final] ≈ 1 atol = testTol
        @test sim.𝒪est[:l2] ≈ 1 atol = testTol
        @test sim.𝒪est[:l∞] ≈ 1 atol = testTol

        sim = test_convergence(
            dts,
            prob,
            QNDF1(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward),
                    function_annotation = Enzyme.Const
                ),
                linsolve = LinearSolve.KrylovJL()
            )
        )
        @test sim.𝒪est[:final] ≈ 1 atol = testTol
        @test sim.𝒪est[:l2] ≈ 1 atol = testTol
        @test sim.𝒪est[:l∞] ≈ 1 atol = testTol
    end

    sim = test_convergence(dts3, prob, QNDF2())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:l2] ≈ 2 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 2 atol = testTol

    sim = test_convergence(dts, prob, QNDF2(nlsolve = NLFunctional()))
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:l2] ≈ 2 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
    @test_nowarn solve(prob, QNDF())

    # MEBDF2
    sim21 = test_convergence(dts, prob, MEBDF2(extrapolant = :linear))
    @test sim21.𝒪est[:final] ≈ 2 atol = testTol

    sim22 = test_convergence(dts, prob, MEBDF2(nlsolve = NLFunctional()), reltol = 1.0e-2)
    @test sim22.𝒪est[:final] ≈ 2 atol = testTol

    sim23 = test_convergence(dts, prob, MEBDF2(nlsolve = NLAnderson()), reltol = 1.0e-2)
    @test sim23.𝒪est[:final] ≈ 2 atol = testTol

    sim24 = test_convergence(
        dts, prob, MEBDF2(nlsolve = NonlinearSolveAlg()), reltol = 1.0e-2
    )
    @test sim24.𝒪est[:final] ≈ 2 atol = testTol

    #FBDF
    @test_nowarn solve(prob, FBDF())
end
