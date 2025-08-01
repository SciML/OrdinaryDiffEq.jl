using OrdinaryDiffEqExplicitRK
using OrdinaryDiffEqCore
using Test
using Printf

"""
Allocation tests for OrdinaryDiffEqExplicitRK solvers.
These tests verify that the step! operation does not allocate during stepping.
"""

@testset "ExplicitRK Allocation Tests" begin
    # Test problem
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    # Based on our testing, these explicit RK solvers are allocation-free
    allocation_free_solvers = [RK4(), BS3(), DP5()]
    
    @testset "Known Allocation-Free Solvers" begin
        for solver in allocation_free_solvers
            @testset "$(typeof(solver)) allocation test" begin
                integrator = init(prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                step!(integrator)  # Setup step may allocate
                
                # Test subsequent steps for zero allocations
                for i in 2:10
                    if integrator.t >= integrator.sol.prob.tspan[2]
                        break
                    end
                    alloc = @allocated step!(integrator)
                    @test alloc == 0
                end
            end
        end
    end
    
    # Test potentially allocating solvers with @test_broken
    @testset "Potentially Allocating Solvers" begin
        # These are marked as broken until allocation issues are resolved
        # (Currently empty as all tested ExplicitRK solvers are allocation-free)
    end
end