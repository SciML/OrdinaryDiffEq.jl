using OrdinaryDiffEqBDF
using OrdinaryDiffEqCore
using Test
using Printf

"""
Allocation tests for OrdinaryDiffEqBDF solvers.
These tests verify that the step! operation should not allocate during stepping.
Currently, many BDF solvers are allocating and marked with @test_broken.
"""

@testset "BDF Allocation Tests" begin
    # Test problem - use a simple linear problem for stiff solvers
    linear_prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    
    # Vector problem
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    vector_prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    # Test known allocating solvers with @test_broken
    allocating_solvers = []
    
    # Try to add available BDF solvers that are currently allocating
    try
        push!(allocating_solvers, QNDF())
    catch
        # QNDF may not be available
    end
    
    try
        push!(allocating_solvers, ABDF2())
    catch
        # ABDF2 may not be available
    end
    
    @testset "Currently Allocating BDF Solvers (@test_broken)" begin
        for solver in allocating_solvers
            @testset "$(typeof(solver)) allocation test (broken)" begin
                integrator = init(linear_prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                step!(integrator)  # Setup step may allocate
                
                # These tests are expected to fail until allocation issues are resolved
                for i in 2:6
                    if integrator.t >= integrator.sol.prob.tspan[2]
                        break
                    end
                    alloc = @allocated step!(integrator)
                    @test_broken alloc == 0
                end
            end
        end
    end
    
    # Placeholder for future allocation-free BDF solvers
    @testset "Future Allocation-Free BDF Solvers" begin
        # When BDF solvers are made allocation-free, move them here from @test_broken
        @test_skip "No allocation-free BDF solvers yet - all currently allocating"
    end
end