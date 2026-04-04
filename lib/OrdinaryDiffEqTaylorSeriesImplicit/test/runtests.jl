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

    sim11 = test_convergence(dts, prob, ImplicitTaylor(μ = 1.0))
    @test sim11.𝒪est[:final]≈1 atol=testTol

    sim12 = test_convergence(dts, prob, ImplicitTaylor(μ = 0.5))
    @test sim12.𝒪est[:final]≈2 atol=testTol

    sim21 = test_convergence(dts, prob, ImplicitTaylor(order = Val(2), μ = 1.0))
    @test sim21.𝒪est[:final]≈2 atol=testTol

    sim22 = test_convergence(dts, prob, ImplicitTaylor(order = Val(2), μ = 0.5))
    @test sim22.𝒪est[:final]≈2 atol=testTol

    # sim23 = test_convergence(dts, prob, ImplicitTaylor(order = Val(2), μ = complex(0.5, √3 / 6)))
    # @test sim23.𝒪est[:final]≈4 atol=testTol

    sim31 = test_convergence(dts, prob, ImplicitTaylor(order = Val(3), μ = 0.5))
    @test sim31.𝒪est[:final]≈4 atol=testTol

    sim41 = test_convergence(dts, prob, ImplicitTaylor(order = Val(4), μ = 0.5))
    @test sim41.𝒪est[:final]≈4 atol=testTol

    # Taylor-Gauss
    simpade1v1 = test_convergence(dts, prob, ImplicitTaylor(order = Val(1), order_q = Val(1)))
    @test simpade1v1.𝒪est[:final]≈2 atol=testTol

    simpade2v2 = test_convergence(dts, prob, ImplicitTaylor(order = Val(2), order_q = Val(2)))
    @test simpade2v2.𝒪est[:final]≈4 atol=testTol

    # Taylor-Radau
    simpade0v1 = test_convergence(dts, prob, ImplicitTaylor(order = Val(0), order_q = Val(1)))
    @test simpade0v1.𝒪est[:final]≈1 atol=testTol

    simpade1v2 = test_convergence(dts, prob, ImplicitTaylor(order = Val(1), order_q = Val(2)))
    @test simpade1v2.𝒪est[:final]≈3 atol=testTol

    # Taylor-Lobatto
    simpade0v2 = test_convergence(dts, prob, ImplicitTaylor(order = Val(0), order_q = Val(2)))
    @test simpade0v2.𝒪est[:final]≈2 atol=testTol

    simpade1v3 = test_convergence(dts, prob, ImplicitTaylor(order = Val(1), order_q = Val(3)))
    @test simpade1v3.𝒪est[:final]≈4 atol=testTol
end
