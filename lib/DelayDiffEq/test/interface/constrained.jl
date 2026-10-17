using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqLowOrderRK
using SciMLBase: ReturnCode
using Test

# Check that numerical solutions approximate analytical solutions,
# independent of problem structure

const alg = MethodOfSteps(BS3(); constrained = true)

# Single constant delay
@testset "single constant delay" begin
    ## Scalar function
    sol_scalar = solve(prob_dde_constant_1delay_scalar, alg; dt = 0.1)

    @test sol_scalar.errors[:l∞] < 3.0e-5
    @test sol_scalar.errors[:final] < 2.1e-5
    @test sol_scalar.errors[:l2] < 1.3e-5

    ## Out-of-place function
    sol_oop = solve(prob_dde_constant_1delay_oop, alg; dt = 0.1)

    @test sol_scalar.t ≈ sol_oop.t && sol_scalar.u ≈ sol_oop[1, :]

    ## In-place function
    sol_ip = solve(prob_dde_constant_1delay_ip, alg; dt = 0.1)

    @test sol_scalar.t ≈ sol_ip.t && sol_scalar.u ≈ sol_ip[1, :]
end

# Two constant delays
@testset "two constant delays" begin
    ## Scalar function
    sol_scalar = solve(prob_dde_constant_2delays_scalar, alg; dt = 0.1)

    @test sol_scalar.errors[:l∞] < 4.1e-6
    @test sol_scalar.errors[:final] < 1.5e-6
    @test sol_scalar.errors[:l2] < 2.3e-6

    ## Out-of-place function
    sol_oop = solve(prob_dde_constant_2delays_oop, alg; dt = 0.1)

    @test sol_scalar.t ≈ sol_oop.t && sol_scalar.u ≈ sol_oop[1, :]

    ## In-place function
    sol_ip = solve(prob_dde_constant_2delays_ip, alg; dt = 0.1)

    @test sol_scalar.t ≈ sol_ip.t && sol_scalar.u ≈ sol_ip[1, :]
end

# No constant lags
@testset "no constant lags" begin
    lag(u, p, t) = 1 + 0.1 * abs(u)
    f(u, h, p, t) = -h(p, t - lag(u, p, t))
    h(p, t) = 1.0
    @testset "$name" for (name, prob) in (
            ("dependent lags", DDEProblem(f, 1.0, h, (0.0, 20.0); dependent_lags = (lag,))),
            ("no lags", DDEProblem(f, 1.0, h, (0.0, 20.0))),
        )
        sol = solve(prob, alg)
        sol_unconstrained = solve(prob, MethodOfSteps(BS3()))

        @test sol.retcode == ReturnCode.Success
        @test sol.t ≈ sol_unconstrained.t && sol.u ≈ sol_unconstrained.u
    end
end
