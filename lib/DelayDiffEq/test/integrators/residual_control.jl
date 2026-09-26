using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqLowOrderRK
using Test

const alg = MethodOfSteps(RK4(); constrained = false)
const prob = prob_dde_constant_1delay_scalar

# reference solution with delays specified
@testset "reference" begin
    sol = solve(prob, alg)

    @test sol.errors[:l∞] < 1.5e-4
    @test sol.errors[:final] < 2.0e-6
    @test sol.errors[:l2] < 5.5e-5
end

# problem without delays specified
const prob_wo = remake(prob; constant_lags = nothing)

# solutions with residual control
@testset "residual control" begin
    sol = solve(prob_wo, alg)

    @test sol.errors[:l∞] < 1.8e-4
    @test sol.errors[:final] < 5.0e-6
    @test sol.errors[:l2] < 9.0e-5

    sol = solve(prob_wo, alg; abstol = 1.0e-9, reltol = 1.0e-6)

    @test sol.errors[:l∞] < 1.0e-7
    @test sol.errors[:final] < 4.5e-9
    @test sol.errors[:l2] < 5.0e-8

    sol = solve(prob_wo, alg; abstol = 1.0e-13, reltol = 1.0e-13)

    # relaxed tests to prevent floating point issues
    @test sol.errors[:l∞] < 5.0e-10
    @test sol.errors[:final] < 4.5e-12
    @test sol.errors[:l2] < 7.7e-11 # 7.7e-12
end

######## Now show that non-residual control is worse
# solutions without residual control
@testset "non-residual control" begin
    sol = solve(prob_wo, MethodOfSteps(OwrenZen5(); constrained = false))

    @test sol.errors[:l∞] > 1.0e-4
    @test sol.errors[:final] > 1.0e-5
    @test sol.errors[:l2] > 5.0e-5

    sol = solve(
        prob_wo, MethodOfSteps(OwrenZen5(); constrained = false); abstol = 1.0e-13,
        reltol = 1.0e-13
    )

    @test sol.errors[:l∞] > 1.0e-4
    @test sol.errors[:final] > 1.0e-5
    @test sol.errors[:l2] > 5.0e-5
end
