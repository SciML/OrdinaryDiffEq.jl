using DelayDiffEq, DDEProblemLibrary
using Test

const alg = MethodOfSteps(BS3())
const prob = prob_dde_constant_1delay_ip

const dde_int = init(prob, alg)
const sol = solve!(dde_int)

@testset "reference" begin
    @test sol.errors[:l∞] < 3.0e-5
    @test sol.errors[:final] < 2.1e-5
    @test sol.errors[:l2] < 1.3e-5
end

@testset "constant lag function" begin
    # constant delay specified as lag function
    prob2 = remake(prob; constant_lags = nothing, dependent_lags = ((u, p, t) -> 1,))
    dde_int2 = init(prob2, alg)
    sol2 = solve!(dde_int2)

    @test getfield.(dde_int.tracked_discontinuities, :t) ≈
        getfield.(dde_int2.tracked_discontinuities, :t)
    @test getfield.(dde_int.tracked_discontinuities, :order) ==
        getfield.(dde_int2.tracked_discontinuities, :order)

    # worse than results above with constant delays specified as scalars
    @test sol2.errors[:l∞] < 3.2e-3
    @test sol2.errors[:final] < 1.2e-4
    @test sol2.errors[:l2] < 1.3e-3

    # simple convergence tests
    @testset "convergence" begin
        sol3 = solve(prob2, alg; abstol = 1.0e-9, reltol = 1.0e-6)

        @test sol3.errors[:l∞] < 3.0e-6
        @test sol3.errors[:final] < 1.4e-7
        @test sol3.errors[:l2] < 9.6e-7

        sol4 = solve(prob2, alg; abstol = 1.0e-13, reltol = 1.0e-13)

        @test sol4.errors[:l∞] < 6.0e-11
        @test sol4.errors[:final] < 4.7e-11
        @test sol4.errors[:l2] < 8.0e-12
    end
end

# without any delays specified is worse
@testset "without delays" begin
    prob2 = remake(prob; constant_lags = nothing)
    sol2 = solve(prob2, alg)

    @test sol2.errors[:l∞] > 5.0e-4
    @test sol2.errors[:final] > 1.0e-6
    @test sol2.errors[:l2] > 2.0e-4
end
