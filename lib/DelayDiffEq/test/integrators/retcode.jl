using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqRosenbrock
using Test
using SciMLBase: ReturnCode

const prob = prob_dde_constant_1delay_ip

@testset for composite in (true, false)
    alg = MethodOfSteps(
        composite ? AutoTsit5(Rosenbrock23()) : Tsit5();
        constrained = false
    )

    sol1 = solve(prob, alg)
    @test sol1.retcode == ReturnCode.Success

    sol2 = solve(prob, alg; maxiters = 1, verbose = false)
    @test sol2.retcode == ReturnCode.MaxIters

    sol3 = solve(prob, alg; dtmin = 5, verbose = false)
    @test sol3.retcode == ReturnCode.DtLessThanMin
end
