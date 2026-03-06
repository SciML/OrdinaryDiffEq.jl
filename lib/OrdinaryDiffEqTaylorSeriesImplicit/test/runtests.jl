using OrdinaryDiffEqTaylorSeriesImplicit, ODEProblemLibrary, DiffEqDevTools, ADTypes
using OrdinaryDiffEqNonlinearSolve: NLFunctional, NLAnderson, NonlinearSolveAlg
using Test, Random
Random.seed!(100)


f_2dlinear = (du, u, p, t) -> (@. du = p * u)
f_2dlinear_analytic = (u0, p, t) -> @. u0 * exp(p * t)
prob_ode_2Dlinear = ODEProblem(
    ODEFunction(f_2dlinear, analytic = f_2dlinear_analytic),
    rand(2), (0.0, 1.0), 1.01
)

## Convergence Testing
testTol = 0.2

@testset "Implicit Solver Convergence Tests" begin # ($(["out-of-place", "in-place"][i]))" for i in 1:2
    # prob = (ODEProblemLibrary.prob_ode_linear,
        # ODEProblemLibrary.prob_ode_2Dlinear)[i]
    prob = ODEProblemLibrary.prob_ode_2Dlinear

    dts = 1 .// 2 .^ (9:-1:5)

    sim11 = test_convergence(dts, prob, ImplicitTaylor(μ = 1.0, extrapolant = :linear))
    @test sim11.𝒪est[:final]≈1 atol=testTol

    sim12 = test_convergence(dts, prob, ImplicitTaylor(μ = 0.5, extrapolant = :linear))
    @test sim12.𝒪est[:final]≈2 atol=testTol

    sim21 = test_convergence(dts, prob, ImplicitTaylor(order = Val(2), μ = 1.0, extrapolant = :linear))
    @test sim21.𝒪est[:final]≈2 atol=testTol

    sim22 = test_convergence(dts, prob, ImplicitTaylor(order = Val(2), μ = 0.5, extrapolant = :linear))
    @test sim22.𝒪est[:final]≈2 atol=testTol

    sim31 = test_convergence(dts, prob, ImplicitTaylor(order = Val(3), μ = 0.5, extrapolant = :linear))
    @test sim31.𝒪est[:final]≈4 atol=testTol

    sim41 = test_convergence(dts, prob, ImplicitTaylor(order = Val(4), μ = 0.5, extrapolant = :linear))
    @test sim41.𝒪est[:final]≈4 atol=testTol
end
