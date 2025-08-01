using OrdinaryDiffEqBDF
using OrdinaryDiffEqCore
using AllocCheck
using Test
using Printf

"""
Allocation tests for OrdinaryDiffEqBDF solvers using AllocCheck.jl.
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
    
    # Test known allocating BDF solvers with @test_broken
    allocating_solvers = [QNDF(), ABDF2()]
    
    @testset "Currently Allocating BDF Solvers (@test_broken)" begin
        for solver in allocating_solvers
            @testset "$(typeof(solver)) allocation check (broken)" begin
                integrator = init(linear_prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                step!(integrator)  # Setup step may allocate
                
                # Use AllocCheck for accurate allocation detection
                allocs = check_allocs(step!, (typeof(integrator),))
                @test_broken length(allocs) == 0  # Should eventually be allocation-free
                
                if length(allocs) > 0
                    println("AllocCheck found $(length(allocs)) allocation sites in $(typeof(solver)) step!:")
                    for (i, alloc) in enumerate(allocs[1:min(3, end)])  # Show first 3
                        println("  $i. $alloc")
                    end
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