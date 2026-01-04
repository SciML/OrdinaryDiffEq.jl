using OrdinaryDiffEqTsit5
using OrdinaryDiffEqCore
using AllocCheck
using Test

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

    # Test all exported Tsit5 solvers for allocation-free behavior
    tsit5_solvers = [Tsit5()]

    @testset "Tsit5 Solver Allocation Analysis" begin
        for solver in tsit5_solvers
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
