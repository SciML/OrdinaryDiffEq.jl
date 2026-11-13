using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqBDF
using OrdinaryDiffEqCore: CompositeAlgorithm
import SciMLBase
using Test

@testset "integrator" begin
    prob = prob_dde_constant_1delay_ip

    integrator = init(prob, MethodOfSteps(AutoTsit5(Rosenbrock23())))
    @test integrator.sol isa SciMLBase.ODESolution
    @test integrator.integrator.sol isa SciMLBase.ODESolution

    sol1 = solve!(integrator)
    @test sol1 isa SciMLBase.ODESolution

    # compare integration grid
    sol2 = solve(prob, MethodOfSteps(Tsit5()))
    @test sol1.t == sol2.t
    @test sol1.u == sol2.u

    # compare interpolation
    @test sol1(0:0.1:1) == sol2(0:0.1:1)
end

@testset "issue #148" begin
    alg = MethodOfSteps(Tsit5())
    compositealg = MethodOfSteps(CompositeAlgorithm((Tsit5(), RK4()), integrator -> 1))

    prob = DDEProblem(
        (du, u, h, p, t) -> (du[1] = -h(p, t - 1)[1]; nothing),
        [1.0], (p, t) -> [1.0], (0.0, 5.0)
    )
    sol = solve(prob, alg)
    compositesol = solve(prob, compositealg)
    @test compositesol isa SciMLBase.ODESolution
    @test sol.t == compositesol.t
    @test sol.u == compositesol.u

    prob_oop = DDEProblem(
        (u, h, p, t) -> -h(p, t - 1), [1.0], (p, t) -> [1.0],
        (0.0, 5.0)
    )
    sol_oop = solve(prob_oop, alg)
    compositesol_oop = solve(prob_oop, compositealg)
    @test compositesol_oop isa SciMLBase.ODESolution
    @test sol_oop.t == compositesol_oop.t
    @test sol_oop.u == compositesol_oop.u
end

@testset "switch to a method with fewer k" begin
    f!(du, u, h, p, t) = (du[1] = -1000 * (u[1] - h(p, t - 1)[1]); nothing)
    f(u, h, p, t) = -1000 .* (u .- h(p, t - 1))
    prob = DDEProblem(f!, [1.0], (p, t) -> [0.5], (0.0, 3.0); constant_lags = [1.0])
    prob_oop = DDEProblem(f, [1.0], (p, t) -> [0.5], (0.0, 3.0); constant_lags = [1.0])
    cb = ContinuousCallback((u, t, integrator) -> t - 2.5, integrator -> (integrator.u[1] += 0.2))

    @testset "Tsit5 to FBDF, iip = $(isinplace(p))" for p in (prob, prob_oop)
        sol = solve(p, MethodOfSteps(AutoTsit5(FBDF())))
        @test sol.retcode == SciMLBase.ReturnCode.Success
    end

    @testset "Tsit5 to Rosenbrock23 with callback, iip = $(isinplace(p))" for p in (prob, prob_oop)
        sol = solve(p, MethodOfSteps(AutoTsit5(Rosenbrock23())); callback = cb)
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test sol.u[end][1] ≈ 0.5 atol = 1.0e-6
    end

    @testset "rejected step after Rosenbrock23 to Tsit5" begin
        alg = CompositeAlgorithm((Rosenbrock23(), Tsit5()), integrator -> integrator.t < 0.3 ? 1 : 2)
        sol = solve(prob, MethodOfSteps(alg))
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test sol.u[end][1] ≈ 0.5 atol = 1.0e-2
    end

    @testset "history extrapolated after a rejected step that switched methods" begin
        g!(du, u, h, p, t) = (du .= -50 .* (u .- h(p, t - 0.01)) .+ 30 * sin(20t); nothing)
        prob_short_lag = DDEProblem(g!, [1.0, 2.0], (p, t) -> [0.5, 0.5], (0.0, 3.0); constant_lags = [0.01])
        alg = CompositeAlgorithm((FBDF(), Tsit5()), integrator -> isodd(integrator.iter) ? 1 : 2)
        sol = solve(prob_short_lag, MethodOfSteps(alg))
        @test sol.retcode == SciMLBase.ReturnCode.Success
    end
end
