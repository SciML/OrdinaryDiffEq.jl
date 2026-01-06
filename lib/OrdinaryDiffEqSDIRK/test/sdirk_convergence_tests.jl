# This definitely needs cleaning
using OrdinaryDiffEqSDIRK, ODEProblemLibrary, DiffEqDevTools, ADTypes
using Test, Random
Random.seed!(100)

## Convergence Testing
testTol = 0.2

@testset "Implicit Solver Convergence Tests ($(["out-of-place", "in-place"][i]))" for i in 1:2
    prob = (
        ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear,
    )[i]

    dts = 1 .// 2 .^ (9:-1:5)

    @show "Very low order"

    sim11 = test_convergence(dts, prob, ImplicitEuler(extrapolant = :linear))
    @test sim11.𝒪est[:final] ≈ 1 atol = testTol

    sim112 = test_convergence(
        dts, prob, ImplicitEuler(nlsolve = NLFunctional()),
        reltol = 1.0e-2
    )
    @test sim112.𝒪est[:final] ≈ 1 atol = testTol

    sim113 = test_convergence(
        dts, prob, ImplicitEuler(nlsolve = NLAnderson()),
        reltol = 1.0e-2
    )
    @test sim113.𝒪est[:final] ≈ 1 atol = testTol

    sim114 = test_convergence(
        dts, prob, ImplicitEuler(nlsolve = NonlinearSolveAlg()),
        reltol = 1.0e-2
    )
    @test sim114.𝒪est[:final] ≈ 1 atol = testTol

    sim13 = test_convergence(dts, prob, ImplicitMidpoint())
    @test sim13.𝒪est[:final] ≈ 2 atol = testTol

    sim132 = test_convergence(dts, prob, ImplicitMidpoint(nlsolve = NLFunctional()))
    @test sim132.𝒪est[:final] ≈ 2 atol = testTol

    sim134 = test_convergence(dts, prob, ImplicitMidpoint(nlsolve = NLAnderson()))
    @test sim134.𝒪est[:final] ≈ 2 atol = testTol

    sim136 = test_convergence(dts, prob, ImplicitMidpoint(nlsolve = NonlinearSolveAlg()))
    @test sim136.𝒪est[:final] ≈ 2 atol = testTol

    sim13 = test_convergence(dts, prob, Trapezoid())
    @test sim13.𝒪est[:final] ≈ 2 atol = testTol

    sim133 = test_convergence(dts, prob, Trapezoid(nlsolve = NLFunctional()))
    @test sim133.𝒪est[:final] ≈ 2 atol = testTol

    sim135 = test_convergence(dts, prob, Trapezoid(nlsolve = NLAnderson()))
    @test sim135.𝒪est[:final] ≈ 2 atol = testTol

    sim137 = test_convergence(dts, prob, Trapezoid(nlsolve = NonlinearSolveAlg()))
    @test sim137.𝒪est[:final] ≈ 2 atol = testTol

    sim14 = test_convergence(dts, prob, TRBDF2())
    @test sim14.𝒪est[:final] ≈ 2 atol = testTol

    sim152 = test_convergence(dts, prob, TRBDF2(autodiff = AutoFiniteDiff()))
    @test sim152.𝒪est[:final] ≈ 2 atol = testTol + 0.1

    sim15 = test_convergence(dts, prob, SDIRK2())
    @test sim15.𝒪est[:final] ≈ 2 atol = testTol

    sim15 = test_convergence(dts, prob, SDIRK22())
    @test sim15.𝒪est[:final] ≈ 2 atol = testTol
    sim152 = test_convergence(dts, prob, SSPSDIRK2())
    @test sim152.𝒪est[:final] ≈ 2 atol = testTol

    @show "Mid SDIRKs"

    sim16 = test_convergence(dts, prob, Kvaerno3())
    @test sim16.𝒪est[:final] ≈ 3 atol = testTol

    sim162 = test_convergence(dts, prob, Kvaerno3(nlsolve = NLFunctional()))
    @test sim162.𝒪est[:final] ≈ 3 atol = testTol

    sim17 = test_convergence(dts, prob, KenCarp3())
    @test sim17.𝒪est[:final] ≈ 3 atol = testTol

    sim019 = test_convergence(dts, prob, CFNLIRK3())
    @test sim019.𝒪est[:final] ≈ 3 atol = testTol

    sim18 = test_convergence(dts, prob, PDIRK44())
    @test sim18.𝒪est[:final] ≈ 4 atol = testTol

    sim182 = test_convergence(dts, prob, PDIRK44(; threading = false))
    @test sim182.𝒪est[:final] ≈ 4 atol = testTol

    dts = (1 / 2) .^ (5:-1:1)
    sim13 = test_convergence(dts, prob, SFSDIRK8())
    @test sim13.𝒪est[:final] ≈ 4 atol = testTol

    dts = (1 / 2) .^ (5:-1:1)
    sim14 = test_convergence(dts, prob, SFSDIRK7())
    @test sim14.𝒪est[:final] ≈ 4 atol = testTol

    dts = (1 / 2) .^ (8:-1:1)
    sim15 = test_convergence(dts, prob, SFSDIRK6())
    @test sim15.𝒪est[:final] ≈ 4 atol = testTol

    dts = (1 / 2) .^ (6:-1:1)
    sim16 = test_convergence(dts, prob, SFSDIRK5())
    @test sim16.𝒪est[:final] ≈ 4 atol = testTol

    dts = 1 .// 2 .^ (7:-1:4)
    sim17 = test_convergence(dts, prob, SFSDIRK4())
    @test sim17.𝒪est[:final] ≈ 4 atol = testTol

    sim18 = test_convergence(dts, prob, Cash4())
    @test sim18.𝒪est[:final] ≈ 4 atol = testTol

    sim19 = test_convergence(dts, prob, Hairer4())
    @test sim19.𝒪est[:final] ≈ 4 atol = testTol

    sim20 = test_convergence(dts, prob, RK46NL())
    @test sim20.𝒪est[:final] ≈ 4 atol = testTol
    @test sim20.𝒪est[:l2] ≈ 4 atol = testTol
    @test sim20.𝒪est[:l∞] ≈ 4 atol = testTol

    sim110 = test_convergence(dts, prob, Hairer42())
    @test sim110.𝒪est[:final] ≈ 4 atol = testTol

    sim111 = test_convergence(dts, prob, Kvaerno4())
    @test sim111.𝒪est[:final] ≈ 4 atol = testTol

    sim112 = test_convergence(dts, prob, KenCarp4())
    @test sim112.𝒪est[:final] ≈ 4 atol = testTol

    sim113 = test_convergence(dts, prob, Kvaerno5())
    @test sim113.𝒪est[:final] ≈ 5 atol = testTol

    sim114 = test_convergence(dts, prob, KenCarp5())
    @test sim114.𝒪est[:final] ≈ 5 atol = testTol

    sim115 = test_convergence(dts, prob, KenCarp5(nlsolve = NLFunctional()))
    @test_broken sim115.𝒪est[:final] ≈ 5 atol = testTol

    sim116 = test_convergence(dts, prob, ESDIRK54I8L2SA())
    @test sim116.𝒪est[:final] ≈ 5 atol = testTol

    sim117 = test_convergence(dts, prob, KenCarp47())
    @test sim117.𝒪est[:final] ≈ 4 atol = testTol

    sim118 = test_convergence(dts, prob, KenCarp58())
    @test sim118.𝒪est[:final] ≈ 5 atol = testTol

    sim119 = test_convergence(dts, prob, ESDIRK436L2SA2())
    @test sim119.𝒪est[:final] ≈ 4 atol = testTol

    sim120 = test_convergence(dts, prob, ESDIRK437L2SA())
    @test sim120.𝒪est[:final] ≈ 4 atol = testTol

    sim121 = test_convergence(dts, prob, ESDIRK547L2SA2())
    @test 5 - testTol <= sim121.𝒪est[:final] <= 6

    dts = (1 / 2) .^ (5:-1:1)
    sim122 = test_convergence(dts, prob, ESDIRK659L2SA())
    @test sim122.𝒪est[:final] ≈ 6 atol = testTol
end
