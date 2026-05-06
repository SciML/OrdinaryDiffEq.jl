using OrdinaryDiffEqTaylorSeriesImplicit, ODEProblemLibrary, DiffEqDevTools, ADTypes
using OrdinaryDiffEqRosenbrock # for reference
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
    # skip linear problems for now, since it can be misleading for stiff solvers
    # prob = ODEProblemLibrary.prob_ode_2Dlinear
    # dts = 1 .// 2 .^ (7:-1:5)

    # sim11 = test_convergence(dts, prob, ImplicitTaylor(μ = 1.0))
    # @test sim11.𝒪est[:final]≈1 atol=testTol

    # sim12 = test_convergence(dts, prob, ImplicitTaylor(μ = 0.5))
    # @test sim12.𝒪est[:final]≈2 atol=testTol

    # sim21 = test_convergence(dts, prob, ImplicitTaylor(order = Val(2), μ = 1.0))
    # @test sim21.𝒪est[:final]≈2 atol=testTol

    # sim22 = test_convergence(dts, prob, ImplicitTaylor(order = Val(2), μ = 0.5))
    # @test sim22.𝒪est[:final]≈2 atol=testTol

    # # sim23 = test_convergence(dts, prob, ImplicitTaylor(order = Val(2), μ = complex(0.5, √3 / 6)))
    # # @test sim23.𝒪est[:final]≈4 atol=testTol

    # sim31 = test_convergence(dts, prob, ImplicitTaylor(order = Val(3), μ = 0.5))
    # @test sim31.𝒪est[:final]≈4 atol=testTol

    # sim41 = test_convergence(dts, prob, ImplicitTaylor(order = Val(4), μ = 0.5))
    # @test sim41.𝒪est[:final]≈4 atol=testTol

    prob_stiff = ODEProblemLibrary.prob_ode_hires
    dts = 2. .^ (-8:-5)
    ref_setup = Dict(:alg => Rodas5P(), :reltol => 1e-14, :abstol => 1e-14)

    # Taylor-Gauss
    sim1v1 = analyticless_test_convergence(dts, prob_stiff, ImplicitTaylor(order = Val(1), order_q = Val(1)), ref_setup)
    @test sim1v1.𝒪est[:final]≈2 atol=testTol

    sim2v2 = analyticless_test_convergence(dts, prob_stiff, ImplicitTaylor(order = Val(2), order_q = Val(2)), ref_setup)
    @test sim2v2.𝒪est[:final]≈4 atol=testTol

    # Taylor-Radau
    sim1v2 = analyticless_test_convergence(dts, prob_stiff, ImplicitTaylor(order = Val(1), order_q = Val(2)), ref_setup)
    @test sim1v2.𝒪est[:final]≈3 atol=testTol

    # Taylor-Lobatto
    sim1v3 = analyticless_test_convergence(dts, prob_stiff, ImplicitTaylor(order = Val(1), order_q = Val(3)), ref_setup)
    @test sim1v3.𝒪est[:final]≈4 atol=testTol
end
