using DelayDiffEq
using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqTsit5
using Test

h(p, t) = 1.0
const tspan = (2.0, 0.0)

f(u, h, p, t) = h(p, t + 1)

# solution on [0,2]
f_analytic(u₀, ::typeof(h), p, t) = t < 1 ? (t^2 - 1) / 2 : t - 1

const dde_f = DDEFunction(f, analytic = f_analytic)

@testset "Without lags" begin
    sol = solve(DDEProblem(dde_f, h, tspan), MethodOfSteps(RK4()))
    @test sol.errors[:l∞] < 1.2e-12 # 1.2e-16
end

@testset "Constant lags" begin
    @testset "incorrect" begin
        prob = DDEProblem(dde_f, h, tspan; constant_lags = [1.0])
        @test_throws ErrorException solve(prob, MethodOfSteps(Tsit5()))
    end

    @testset "correct" begin
        prob = DDEProblem(dde_f, h, tspan; constant_lags = [-1.0])

        dde_int = init(prob, MethodOfSteps(Tsit5()))
        @test dde_int.tracked_discontinuities == [Discontinuity(-2.0, 1)]
        @test dde_int.d_discontinuities_propagated.valtree == [Discontinuity(-1.0, 2)]
        @test isempty(dde_int.opts.d_discontinuities)

        sol = solve!(dde_int)
        @test sol.errors[:l∞] < 3.9e-12 # 3.9e-15
        @test dde_int.tracked_discontinuities ==
            [Discontinuity(-2.0, 1), Discontinuity(-1.0, 2)]
        @test isempty(dde_int.d_discontinuities_propagated)
        @test isempty(dde_int.opts.d_discontinuities)
    end
end

@testset "Dependent lags" begin
    prob = DDEProblem(dde_f, h, tspan; dependent_lags = ((u, p, t) -> -1.0,))

    dde_int = init(prob, MethodOfSteps(Tsit5()))
    @test isempty(dde_int.opts.d_discontinuities)

    sol = solve!(dde_int)
    @test sol.errors[:l∞] < 3.9e-12 # 3.9e-15
    @test dde_int.tracked_discontinuities == [Discontinuity(-2.0, 1)]
end

@testset "dt and dtmax" begin
    prob = DDEProblem(dde_f, h, tspan)

    dde_int = init(prob, MethodOfSteps(RK4()); dt = 0.1, dtmax = 0.5)
    @test dde_int.dt == -0.1
    @test dde_int.opts.dtmax == -0.5

    sol1 = solve!(dde_int)
    sol2 = solve(prob, MethodOfSteps(RK4()); dt = -0.1, dtmax = -0.5)
    @test sol1.t == sol2.t
    @test sol1.u == sol2.u
end
