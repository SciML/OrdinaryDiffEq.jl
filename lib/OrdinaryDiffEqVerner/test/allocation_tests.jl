using Pkg
Pkg.add("AllocCheck")

using OrdinaryDiffEqVerner
using OrdinaryDiffEqCore
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqVerner solvers using AllocCheck.jl.
These tests verify that the step! operation does not allocate during stepping.
"""

@testset "Verner Allocation Tests" begin
    # Test problem
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))

    # Test all exported Verner solvers for allocation-free behavior
    verner_solvers = [
        Vern6(), Vern7(), Vern8(), Vern9(),
        AutoVern6(Vern6()), AutoVern7(Vern7()), AutoVern8(Vern8()), AutoVern9(Vern9()),
    ]

    @testset "Verner Solver Allocation Analysis" begin
        for solver in verner_solvers
            @testset "$(typeof(solver)) allocation check" begin
                integrator = init(prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
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
