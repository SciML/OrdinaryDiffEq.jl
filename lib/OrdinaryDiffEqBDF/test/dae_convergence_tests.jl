using OrdinaryDiffEqBDF, DiffEqDevTools, ADTypes
using Test, Random
Random.seed!(100)

dts = 1 .// 2 .^ (9:-1:5)
testTol = 0.2

f_dae_linear = (res, du, u, p, t) -> (@. res = du - u)
function f_dae_linear_jac(J, du, u, p, gamma, t)
    J[1, 1] = gamma - 1.0
    return J[2, 2] = gamma - 1.0
end
f_dae_linear_analytic = (du0, u0, p, t) -> @. u0 * exp(t)
prob_dae_linear_iip = DAEProblem(
    DAEFunction(
        f_dae_linear;
        analytic = f_dae_linear_analytic
    ),
    [1.0, 1.0], [1.0, 1.0], (0.0, 1.0)
)

@testset "DAE Solver Convergence Tests (in-place, ad jac)" begin
    prob = prob_dae_linear_iip

    sim11 = test_convergence(dts, prob, DImplicitEuler())
    @test sim11.𝒪est[:final] ≈ 1 atol = testTol

    sim12 = test_convergence(dts, prob, DImplicitEuler(; autodiff = AutoFiniteDiff()))
    @test sim12.𝒪est[:final] ≈ 1 atol = testTol

    sim13 = test_convergence(dts, prob, DABDF2())
    @test sim13.𝒪est[:final] ≈ 2 atol = testTol

    sim14 = test_convergence(dts, prob, DABDF2(; autodiff = AutoFiniteDiff()))
    @test sim14.𝒪est[:final] ≈ 2 atol = testTol

    @test_nowarn solve(prob, DFBDF())
end

prob_dae_linear_iip_jac = DAEProblem(
    DAEFunction(
        f_dae_linear;
        jac = f_dae_linear_jac,
        analytic = f_dae_linear_analytic
    ),
    [1.0, 1.0], [1.0, 1.0], (0.0, 1.0)
)

@testset "DAE Solver Convergence Tests (in-place, custom jac)" begin
    prob = prob_dae_linear_iip

    sim11 = test_convergence(dts, prob, DImplicitEuler())
    @test sim11.𝒪est[:final] ≈ 1 atol = testTol

    sim12 = test_convergence(dts, prob, DImplicitEuler(; autodiff = AutoFiniteDiff()))
    @test sim12.𝒪est[:final] ≈ 1 atol = testTol

    sim13 = test_convergence(dts, prob, DABDF2())
    @test sim13.𝒪est[:final] ≈ 2 atol = testTol

    sim14 = test_convergence(dts, prob, DABDF2(; autodiff = AutoFiniteDiff()))
    @test sim14.𝒪est[:final] ≈ 2 atol = testTol

    @test_nowarn solve(prob, DFBDF())
end

# IIP separated Jacobians: jac_u and jac_du
function f_dae_linear_jac_u_iip(J, du, u, p, t)
    fill!(J, 0)
    J[1, 1] = -1.0
    J[2, 2] = -1.0
    nothing
end
function f_dae_linear_jac_du_iip(J, du, u, p, t)
    fill!(J, 0)
    J[1, 1] = 1.0
    J[2, 2] = 1.0
    nothing
end
prob_dae_linear_iip_sep = DAEProblem(
    DAEFunction(
        (res, du, u, p, t) -> (@. res = du - u);
        jac_u = f_dae_linear_jac_u_iip,
        jac_du = f_dae_linear_jac_du_iip,
        analytic = (du0, u0, p, t) -> @. u0 * exp(t)
    ),
    [1.0, 1.0], [1.0, 1.0], (0.0, 1.0)
)

@testset "DAE Solver Convergence Tests (in-place, separated jac_u/jac_du)" begin
    prob = prob_dae_linear_iip_sep

    sim11 = test_convergence(dts, prob, DImplicitEuler())
    @test sim11.𝒪est[:final] ≈ 1 atol = testTol

    sim13 = test_convergence(dts, prob, DABDF2())
    @test sim13.𝒪est[:final] ≈ 2 atol = testTol

    @test_nowarn solve(prob, DFBDF())
end

f_dae_linear = (du, u, p, t) -> (@. du - u)
function f_dae_linear_jac(du, u, p, gamma, t)
    J = zeros(2, 2)
    J[1, 1] = gamma - 1.0
    return J[2, 2] = gamma - 1.0
end
f_dae_linear_analytic = (du0, u0, p, t) -> @. u0 * exp(t)
prob_dae_linear_oop = DAEProblem(
    DAEFunction(
        f_dae_linear;
        analytic = f_dae_linear_analytic
    ),
    1.0, 1.0, (0.0, 1.0)
)

@testset "DAE Solver Convergence Tests (out-of-place, ad jac)" begin
    prob = prob_dae_linear_oop

    sim21 = test_convergence(dts, prob, DImplicitEuler())
    @test sim21.𝒪est[:final] ≈ 1 atol = testTol

    sim22 = test_convergence(dts, prob, DImplicitEuler(; autodiff = AutoFiniteDiff()))
    @test sim22.𝒪est[:final] ≈ 1 atol = testTol

    sim23 = test_convergence(dts, prob, DABDF2())
    @test sim23.𝒪est[:final] ≈ 2 atol = testTol

    sim24 = test_convergence(dts, prob, DABDF2(; autodiff = AutoFiniteDiff()))
    @test sim24.𝒪est[:final] ≈ 2 atol = testTol

    @test_nowarn solve(prob, DFBDF())
end

prob_dae_linear_oop = DAEProblem(
    DAEFunction(
        f_dae_linear;
        jac = f_dae_linear_jac,
        analytic = f_dae_linear_analytic
    ),
    1.0, 1.0, (0.0, 1.0)
)

@testset "DAE Solver Convergence Tests (out-of-place, custom jac)" begin
    prob = prob_dae_linear_oop

    sim21 = test_convergence(dts, prob, DImplicitEuler())
    @test sim21.𝒪est[:final] ≈ 1 atol = testTol

    sim22 = test_convergence(dts, prob, DImplicitEuler(; autodiff = AutoFiniteDiff()))
    @test sim22.𝒪est[:final] ≈ 1 atol = testTol

    sim23 = test_convergence(dts, prob, DABDF2())
    @test sim23.𝒪est[:final] ≈ 2 atol = testTol

    sim24 = test_convergence(dts, prob, DABDF2(; autodiff = AutoFiniteDiff()))
    @test sim24.𝒪est[:final] ≈ 2 atol = testTol

    @test_nowarn solve(prob, DFBDF())
end

# OOP separated Jacobians: jac_u and jac_du (scalar problem)
f_dae_linear_jac_u_oop(du, u, p, t) = -1.0
f_dae_linear_jac_du_oop(du, u, p, t) = 1.0
prob_dae_linear_oop_sep = DAEProblem(
    DAEFunction(
        (du, u, p, t) -> (@. du - u);
        jac_u = f_dae_linear_jac_u_oop,
        jac_du = f_dae_linear_jac_du_oop,
        analytic = (du0, u0, p, t) -> @. u0 * exp(t)
    ),
    1.0, 1.0, (0.0, 1.0)
)

@testset "DAE Solver Convergence Tests (out-of-place, separated jac_u/jac_du)" begin
    prob = prob_dae_linear_oop_sep

    sim21 = test_convergence(dts, prob, DImplicitEuler())
    @test sim21.𝒪est[:final] ≈ 1 atol = testTol

    sim23 = test_convergence(dts, prob, DABDF2())
    @test sim23.𝒪est[:final] ≈ 2 atol = testTol

    @test_nowarn solve(prob, DFBDF())
end
