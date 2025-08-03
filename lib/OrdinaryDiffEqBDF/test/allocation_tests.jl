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
    # Test problem - use a simple linear problem for stiff solvers
    linear_prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    
    # Vector problem
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    vector_prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    # Test all exported BDF solvers for allocation-free behavior
    bdf_solvers = [ABDF2(), QNDF1(), QBDF1(), QNDF2(), QBDF2(), QNDF(), QBDF(), FBDF(),
                   SBDF(), SBDF2(), SBDF3(), SBDF4(), MEBDF2(), IMEXEuler(), IMEXEulerARK(),
                   DABDF2(), DImplicitEuler(), DFBDF()]
    
    @testset "BDF Solver Allocation Analysis" begin
        for solver in bdf_solvers
            @testset "$(typeof(solver)) allocation check" begin
                integrator = init(linear_prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                step!(integrator)  # Setup step may allocate
                
                # Use AllocCheck for accurate allocation detection
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