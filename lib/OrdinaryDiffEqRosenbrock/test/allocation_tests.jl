using Pkg
Pkg.add("AllocCheck")

using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqCore
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqRosenbrock solvers using AllocCheck.jl.
These tests verify that the step! operation should not allocate during stepping.
Currently, Rosenbrock solvers are allocating and marked with @test_broken.
"""

@testset "Rosenbrock Allocation Tests" begin
    # Test problem - use a simple linear problem for stiff solvers
    linear_prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))

    # Vector problem
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    vector_prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))

    # Test all exported Rosenbrock solvers for allocation-free behavior
    rosenbrock_solvers = [
        Rosenbrock23(), Rosenbrock32(), RosShamp4(), Veldd4(), Velds4(), GRK4T(), GRK4A(),
        Rodas3(), Rodas23W(), Rodas3P(), Rodas4(), Rodas42(), Rodas4P(), Rodas4P2(), Rodas5(),
        Rodas5P(), Rodas5Pe(), Rodas5Pr(), Rodas6P(),
    ]

    @testset "Rosenbrock Solver Allocation Analysis" begin
        for solver in rosenbrock_solvers
            @testset "$(typeof(solver)) allocation check" begin
                integrator = init(
                    linear_prob, solver, dt = 0.1,
                    save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6
                )
                step!(integrator)  # Setup step may allocate

                # Use AllocCheck for accurate allocation detection
                allocs = check_allocs(step!, (typeof(integrator),))

                # These solvers should be allocation-free, but mark as broken for now
                # to verify with AllocCheck (more accurate than @allocated)
                @test_broken length(allocs) == 0

                # However, we don't want to allow any dynamic dispatch from this module
                @test count(
                    t -> t isa AllocCheck.DynamicDispatch &&
                        any(s -> contains(string(s.file), "Rosenbrock"), t.backtrace),
                    allocs
                ) == 0

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
