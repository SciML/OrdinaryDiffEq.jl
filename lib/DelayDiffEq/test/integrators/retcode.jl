using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqCore: DEVerbosity
using SciMLLogging: None
using Test
using SciMLBase: ReturnCode

const prob = prob_dde_constant_1delay_ip

@testset "composite: $composite" for composite in (true, false)
    alg = MethodOfSteps(
        composite ? AutoTsit5(Rosenbrock23()) : Tsit5();
        constrained = false
    )

    sol1 = solve(prob, alg)
    @test sol1.retcode == ReturnCode.Success

    sol2 = solve(prob, alg; maxiters = 1, verbose = DEVerbosity(None()))
    @test sol2.retcode == ReturnCode.MaxIters

    sol3 = solve(prob, alg; dtmin = 5, verbose = DEVerbosity(None()))
    @test sol3.retcode == ReturnCode.DtLessThanMin
end
