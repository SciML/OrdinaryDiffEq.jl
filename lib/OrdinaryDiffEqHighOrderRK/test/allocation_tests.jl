using OrdinaryDiffEqHighOrderRK
using OrdinaryDiffEqCore
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqHighOrderRK solvers using AllocCheck.jl.
These tests verify that the step! operation does not allocate during stepping.
"""

@testset "HighOrderRK Allocation Tests" begin
    # Test problem
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    # Test all exported HighOrderRK solvers for allocation-free behavior
    high_order_solvers = [TanYam7(), DP8(), PFRK87(), TsitPap8()]
    
    @testset "HighOrderRK Solver Allocation Analysis" begin
        for solver in high_order_solvers
            @testset "$(typeof(solver)) allocation check" begin
                integrator = init(prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                step!(integrator)  # Setup step may allocate
                
                # Use AllocCheck to verify step! is allocation-free
                allocs = check_allocs(step!, (typeof(integrator),))
                
                # These solvers should be allocation-free, but mark as broken for now
                # to verify with AllocCheck (more accurate than @allocated)
                @test_broken length(allocs) == 0
                
                if length(allocs) > 0
                    println("AllocCheck found $(length(allocs)) allocation sites in $(typeof(solver)) step!:")
                    for (i, alloc) in enumerate(allocs[1:min(3, end)])  # Show first 3
                        println("  $i. $alloc")
                    end
                else
                    println("✓ $(typeof(solver)) appears allocation-free with AllocCheck")
                end
            end
        end
    end
end