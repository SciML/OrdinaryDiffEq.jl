using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqTsit5
using Test

const prob = prob_dde_constant_1delay_long_ip

@testset for constrained in (false, true)
    sol = solve(prob, MethodOfSteps(Tsit5(); constrained = constrained))

    @test allunique(sol.t)
end
