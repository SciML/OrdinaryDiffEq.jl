using OrdinaryDiffEqTsit5
using OrdinaryDiffEqCore
using AllocCheck
using Test
using Printf

"""
Allocation tests for OrdinaryDiffEqTsit5 solvers using AllocCheck.jl.
These tests verify that the step! operation does not allocate during stepping.
"""

@testset "Tsit5 Allocation Tests" begin
    # Test problem
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    @testset "Tsit5 Step Allocation Analysis" begin
        integrator = init(prob, Tsit5(), save_everystep=false, abstol=1e-6, reltol=1e-6)
        
        # First step may allocate for setup
        step!(integrator)
        
        # Use AllocCheck to verify step! is allocation-free
        allocs = check_allocs(step!, (typeof(integrator),))
        
        # Currently expect allocations (mark as broken until fixed)
        @test_broken length(allocs) == 0
        
        if length(allocs) > 0
            println("AllocCheck found $(length(allocs)) allocation sites in Tsit5 step!:")
            for (i, alloc) in enumerate(allocs[1:min(3, end)])  # Show first 3
                println("  $i. $alloc")
            end
        end
    end
    
    @testset "Current Solver Status" begin
        # Test currently allocating solvers with @test_broken until they're fixed
        allocating_solvers = [Tsit5()]
        
        for solver in allocating_solvers
            @testset "$(typeof(solver)) allocation check (broken)" begin
                integrator = init(prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                step!(integrator)  # Setup step
                
                # Use AllocCheck for accurate allocation detection
                allocs = check_allocs(step!, (typeof(integrator),))
                @test_broken length(allocs) == 0  # Should eventually be allocation-free
            end
        end
    end
end