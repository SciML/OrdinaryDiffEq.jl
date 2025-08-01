using OrdinaryDiffEqTsit5
using OrdinaryDiffEqCore
using Test
using Printf

"""
Allocation tests for OrdinaryDiffEqTsit5 solvers.
These tests verify that the step! operation does not allocate during stepping.
"""

@testset "Tsit5 Allocation Tests" begin
    # Test problem
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    @testset "Tsit5 Step Allocations (@test_broken)" begin
        integrator = init(prob, Tsit5(), save_everystep=false, abstol=1e-6, reltol=1e-6)
        
        # First step may allocate for setup
        step!(integrator)
        
        # Subsequent steps should eventually be allocation-free (currently broken)
        for i in 2:10
            if integrator.t >= integrator.sol.prob.tspan[2]
                break
            end
            
            alloc = @allocated step!(integrator)
            @test_broken alloc == 0  # Currently allocating - needs to be fixed
        end
    end
    
    # Test currently allocating solvers with @test_broken until they're fixed
    @testset "Currently Allocating Solvers (@test_broken)" begin
        allocating_solvers = [Tsit5()]
        
        for solver in allocating_solvers
            @testset "$(typeof(solver)) allocation test (broken)" begin
                integrator = init(prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                step!(integrator)  # Setup step
                
                # Test 5 subsequent steps - these should eventually be allocation-free
                for i in 2:6
                    if integrator.t >= integrator.sol.prob.tspan[2]
                        break
                    end
                    alloc = @allocated step!(integrator)
                    @test_broken alloc == 0  # Mark as broken until allocation issues are resolved
                end
            end
        end
    end
end