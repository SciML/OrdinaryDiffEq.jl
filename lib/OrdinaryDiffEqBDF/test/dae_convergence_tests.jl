using OrdinaryDiffEqBDF, DiffEqDevTools, ADTypes
using SciMLBase: successful_retcode
import SciMLBase
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
    return nothing
end
function f_dae_linear_jac_du_iip(J, du, u, p, t)
    fill!(J, 0)
    J[1, 1] = 1.0
    J[2, 2] = 1.0
    return nothing
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

# Regression test for issue #3960: DABDF2 adaptive error estimate blew up on the
# first BDF2 step (the FSAL startup left integrator.fsallast stale, so fsalfirst
# was ~0 at the first estimate), forcing dt below eps and returning Unstable on a
# trivial index-1 DAE that DImplicitEuler/DFBDF solve. The fixed-dt convergence
# tests above do not exercise the adaptive controller, so this uses adaptive solves.
@testset "DABDF2 adaptive index-1 DAE (issue #3960)" begin
    # u1' = -u1 ;  0 = u2 - u1     exact: u1 = u2 = exp(-t)
    f_index1_iip = (res, du, u, p, t) -> begin
        res[1] = du[1] + u[1]
        res[2] = u[2] - u[1]
        return nothing
    end
    prob_iip = DAEProblem(
        f_index1_iip, [-1.0, -1.0], [1.0, 1.0], (0.0, 1.0);
        differential_vars = [true, false]
    )
    # scalar out-of-place hits the DABDF2ConstantCache path
    f_scalar_oop = (du, u, p, t) -> du + u
    prob_oop = DAEProblem(f_scalar_oop, -1.0, 1.0, (0.0, 1.0))

    exact = exp(-1.0)
    prev_err = Inf
    for tol in (1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10)
        # atol 2e-3 cleanly separates a correct solve (worst observed err ~8e-4)
        # from the pre-fix Unstable behavior (u stuck near the u0 = 1.0, err ~0.6)
        sol_iip = solve(prob_iip, DABDF2(), abstol = tol, reltol = tol)
        @test SciMLBase.successful_retcode(sol_iip)
        @test sol_iip.u[end][1] ≈ exact atol = 2.0e-3
        @test sol_iip.u[end][2] ≈ exact atol = 2.0e-3

        sol_oop = solve(prob_oop, DABDF2(), abstol = tol, reltol = tol)
        @test SciMLBase.successful_retcode(sol_oop)
        @test sol_oop.u[end] ≈ exact atol = 2.0e-3

        # tightening the tolerance must improve accuracy (estimator is meaningful)
        err = abs(sol_iip.u[end][1] - exact)
        @test err <= prev_err
        prev_err = err
    end
end
