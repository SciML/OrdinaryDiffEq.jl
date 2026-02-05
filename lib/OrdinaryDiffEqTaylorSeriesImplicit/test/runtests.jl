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

@testset "Implicit Solver Convergence Tests ($(["out-of-place", "in-place"][i]))" for i in 1:2
    prob = (ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear)[i]

    dts = 1 .// 2 .^ (9:-1:5)

    μ = 1.0

    @show "Very low order"

    sim11 = test_convergence(dts, prob, ImplicitTaylor1(μ = μ, extrapolant = :linear))
    @test sim11.𝒪est[:final]≈1 atol=testTol

    sim112 = test_convergence(dts, prob, ImplicitTaylor1(μ = μ, nlsolve = NLFunctional()),
        reltol = 1e-2)
    @test sim112.𝒪est[:final]≈1 atol=testTol

    sim113 = test_convergence(dts, prob, ImplicitTaylor1(μ = μ, nlsolve = NLAnderson()),
        reltol = 1e-2)
    @test sim113.𝒪est[:final]≈1 atol=testTol

    sim114 = test_convergence(dts, prob, ImplicitTaylor1(μ = μ, nlsolve = NonlinearSolveAlg()),
        reltol = 1e-2)
    @test sim114.𝒪est[:final]≈1 atol=testTol
end
