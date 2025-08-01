using OrdinaryDiffEqHighOrderRK
using OrdinaryDiffEqCore
using Test
using Printf

"""
Allocation tests for OrdinaryDiffEqHighOrderRK solvers.
These tests verify that the step! operation does not allocate during stepping.
"""

@testset "HighOrderRK Allocation Tests" begin
    # Test problem
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    # Based on our testing, these high-order RK solvers are allocation-free
    allocation_free_solvers = [Vern6(), Vern7(), Vern8(), Vern9()]
    
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
end