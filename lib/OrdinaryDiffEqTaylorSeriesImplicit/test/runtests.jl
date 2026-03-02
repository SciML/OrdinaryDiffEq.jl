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

    μ = 1.0

    # sol = solve(prob, ImplicitTaylor1(μ = μ, extrapolant = :linear))

    # @show sol.stats

    sim11 = test_convergence(dts, prob, ImplicitTaylor1(μ = μ, extrapolant = :linear))
    @test sim11.𝒪est[:final]≈1 atol=testTol
end
