using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqDefault
using OrdinaryDiffEqCore: CompositeAlgorithm
using OrdinaryDiffEqNonlinearSolve: NLAnderson
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

@testset "AutoSwitch keeps constrained and fpsolve" begin
    # u(t) = exp(t) solves this problem, and its unconstrained steps are longer than the lag
    τ = 0.1
    prob_exp = DDEProblem(
        (u, h, p, t) -> exp(τ) * h(p, t - τ), 1.0, (p, t) -> exp(t), (0.0, 2.0);
        constant_lags = [τ]
    )
    prob_overlap = remake(
        prob_exp; constant_lags = [0.25], dependent_lags = ((u, p, t) -> τ,)
    )

    for alg in (AutoTsit5(Rosenbrock23()), DefaultODEAlgorithm())
        fpsolve = NLAnderson()
        integrator = init(
            prob_overlap, MethodOfSteps(alg; constrained = true, fpsolve);
            adaptive = false, dt = 0.25, abstol = 1.0e-9, reltol = 1.0e-9
        )
        @test integrator.opts.dtmax == 0.25
        @test integrator.fpsolver.alg === fpsolve
        sol = solve!(integrator)
        @test SciMLBase.successful_retcode(sol)
        @test sol.stats.nsolve > 0
        @test maximum(abs(sol(t) - exp(t)) for t in 0.0:0.125:2.0) < 1.0e-6

        sol = solve(prob_exp, MethodOfSteps(alg; constrained = true))
        @test maximum(diff(sol.t)) ≈ τ
        @test sol.u ≈ exp.(sol.t) rtol = 1.0e-8

        sol = solve(prob_exp, MethodOfSteps(alg; constrained = true); adaptive = false, dt = 0.3)
        @test SciMLBase.successful_retcode(sol)
        @test maximum(diff(sol.t)) ≈ τ
        @test sol.stats.nfpiter == 0
    end
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
