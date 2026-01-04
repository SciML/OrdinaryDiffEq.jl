using OrdinaryDiffEqAdamsBashforthMoulton, ODEProblemLibrary
import OrdinaryDiffEqCore: OrdinaryDiffEqAdaptiveAlgorithm
using Test
using Static

algorithms = [
    AB3, AB4, AB5, ABM32, ABM43, ABM54, VCAB3, VCAB4, VCAB5, VCABM3, VCABM4, VCABM5, VCABM,
]

problem = ODEProblemLibrary.prob_ode_linear

@testset "Regression test for threading versions vs non threading versions" begin
    @testset "$ALG" for ALG in algorithms
        if ALG isa OrdinaryDiffEqAdaptiveAlgorithm
            sol_thread = solve(problem, ALG(Static.True()))
            sol_nothread = solve(problem, ALG(Static.False()))
        else
            sol_thread = solve(problem, ALG(Static.True()), dt = 1 // 2^9)
            sol_nothread = solve(problem, ALG(Static.False()), dt = 1 // 2^9)
        end
        @test all(sol_nothread .== sol_thread)
    end
end
