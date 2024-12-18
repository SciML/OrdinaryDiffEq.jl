using OrdinaryDiffEqAdamsBashforthMoulton, ODEProblemLibrary
using Test
using Static

free_timestep_algorithms = [ VCAB3, VCAB4, VCAB5, VCABM3, VCABM4, VCABM5, VCABM]

fixed_timestep_algorithms = [AB3,AB4, AB5, ABM32, ABM43, ABM54]

problem = ODEProblemLibrary.prob_ode_linear

function test_alg(ALG; kwargs...)
   sol_thread = solve(problem, ALG(Static.True()); kwargs...)
   sol_nothread = solve(problem, ALG(Static.False()); kwargs...)

   @test all(sol_nothread .== sol_thread)
end


@testset "Regression test for threading versions vs non threading versions" begin
    @testset "$ALG" for ALG in fixed_timestep_algorithms
        test_alg(ALG, dt=1//2^9)
    end
    @testset "$ALG" for ALG in free_timestep_algorithms
        test_alg(ALG)
    end

end
