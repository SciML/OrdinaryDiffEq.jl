using OrdinaryDiffEq, DiffEqDevTools
using Test, Random
Random.seed!(100)

dts = 1 .// 2 .^ (9:-1:5)
testTol = 0.2

f_dae_linear = (res, du, u, p, t) -> (@. res = du - u)
f_dae_linear_analytic = (du0, u0, p, t) -> @. u0 * exp(t)
prob_dae_linear_iip = DAEProblem(
    DAEFunction(f_dae_linear;
        analytic = f_dae_linear_analytic),
    [1.0, 1.0], [1.0, 1.0], (0.0, 1.0))

@testset "DAE Solver Convergence Tests (in-place)" begin
    prob = prob_dae_linear_iip

    sim11 = test_convergence(dts, prob, DImplicitEuler())
    @test sim11.𝒪est[:final]≈1 atol=testTol

    sim12 = test_convergence(dts, prob, DImplicitEuler(; autodiff = false))
    @test sim12.𝒪est[:final]≈1 atol=testTol

    sim13 = test_convergence(dts, prob, DABDF2())
    @test sim13.𝒪est[:final]≈2 atol=testTol

    sim14 = test_convergence(dts, prob, DABDF2(; autodiff = false))
    @test sim14.𝒪est[:final]≈2 atol=testTol

    @test_nowarn solve(prob, DFBDF())
end

f_dae_linear = (du, u, p, t) -> (@. du - u)
f_dae_linear_analytic = (du0, u0, p, t) -> @. u0 * exp(t)
prob_dae_linear_oop = DAEProblem(
    DAEFunction(f_dae_linear;
        analytic = f_dae_linear_analytic),
    1.0, 1.0, (0.0, 1.0))

@testset "DAE Solver Convergence Tests (out-of-place)" begin
    prob = prob_dae_linear_oop

    sim21 = test_convergence(dts, prob, DImplicitEuler())
    @test sim21.𝒪est[:final]≈1 atol=testTol

    sim22 = test_convergence(dts, prob, DImplicitEuler(; autodiff = false))
    @test sim22.𝒪est[:final]≈1 atol=testTol

    sim23 = test_convergence(dts, prob, DABDF2())
    @test sim23.𝒪est[:final]≈2 atol=testTol

    sim24 = test_convergence(dts, prob, DABDF2(; autodiff = false))
    @test sim24.𝒪est[:final]≈2 atol=testTol

    @test_nowarn solve(prob, DFBDF())
end