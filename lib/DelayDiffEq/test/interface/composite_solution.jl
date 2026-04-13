using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqLowOrderRK
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
