using OrdinaryDiffEqBDF
using OrdinaryDiffEqCore
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqBDF solvers using AllocCheck.jl.
These tests verify that the step! operation should not allocate during stepping.
Currently, many BDF solvers are allocating and marked with @test_broken.
"""

@testset "BDF Allocation Tests" begin
    # Test problem - use a simple linear problem for standard BDF solvers
    linear_prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    
    # Vector problem for in-place methods
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    vector_prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    # Split problem for IMEX methods
    function f1!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = 0.0
    end
    function f2!(du, u, p, t)
        du[1] = 0.0
        du[2] = -1.5 * u[2]
    end
    split_prob = SplitODEProblem(f1!, f2!, [1.0, 1.0], (0.0, 1.0))
    
    # Group 1: Standard BDF methods (use linear problem)
    standard_bdf_solvers = [ABDF2(), QNDF1(), QBDF1(), QNDF2(), QBDF2(), QNDF(), QBDF(), FBDF(), MEBDF2()]
    
    # Group 2: SBDF methods (use linear problem)
    sbdf_solvers = [SBDF(order=2), SBDF2(), SBDF3(), SBDF4()]
    
    # Group 3: IMEX methods (use split problem)
    imex_solvers = [IMEXEuler(), IMEXEulerARK()]
    
    # Group 4: Dual/DAE methods (use vector problem)
    dual_solvers = [DABDF2(), DImplicitEuler(), DFBDF()]
    
    @testset "Standard BDF Solver Allocation Analysis" begin
        for solver in standard_bdf_solvers
            @testset "$(typeof(solver)) allocation check" begin
                integrator = init(linear_prob, solver, dt=0.1, save_everystep=false, abstol=1e-6, reltol=1e-6)
                step!(integrator)  # Setup step may allocate
                
                # Use AllocCheck for accurate allocation detection
                allocs = check_allocs(step!, (typeof(integrator),))
                
                # These solvers should be allocation-free, but mark as broken for now
                @test length(allocs) == 0 broken=true
                
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
    
    @testset "SBDF Solver Allocation Analysis" begin
        for solver in sbdf_solvers
            @testset "$(typeof(solver)) allocation check" begin
                integrator = init(linear_prob, solver, dt=0.1, save_everystep=false, abstol=1e-6, reltol=1e-6)
                step!(integrator)  # Setup step may allocate
                
                # Use AllocCheck for accurate allocation detection
                allocs = check_allocs(step!, (typeof(integrator),))
                
                # These solvers should be allocation-free, but mark as broken for now
                @test length(allocs) == 0 broken=true
                
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
    
    @testset "IMEX Solver Allocation Analysis" begin
        for solver in imex_solvers
            @testset "$(typeof(solver)) allocation check" begin
                integrator = init(split_prob, solver, dt=0.1, save_everystep=false, abstol=1e-6, reltol=1e-6)
                step!(integrator)  # Setup step may allocate
                
                # Use AllocCheck for accurate allocation detection
                allocs = check_allocs(step!, (typeof(integrator),))
                
                # These solvers should be allocation-free, but mark as broken for now
                @test length(allocs) == 0 broken=true
                
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
    
    @testset "Dual/DAE Solver Allocation Analysis" begin
        for solver in dual_solvers
            @testset "$(typeof(solver)) allocation check" begin
                integrator = init(vector_prob, solver, dt=0.1, save_everystep=false, abstol=1e-6, reltol=1e-6)
                step!(integrator)  # Setup step may allocate
                
                # Use AllocCheck for accurate allocation detection
                allocs = check_allocs(step!, (typeof(integrator),))
                
                # These solvers should be allocation-free, but mark as broken for now
                @test length(allocs) == 0 broken=true
                
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