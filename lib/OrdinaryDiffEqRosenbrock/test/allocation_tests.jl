using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqCore
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqRosenbrock solvers using AllocCheck.jl.
These tests verify that perform_step! (the core stepping function) should not allocate.
Note: We test perform_step! directly rather than step! because step! includes saving
operations that will naturally allocate. The core numerical stepping should be allocation-free.
"""

@testset "Rosenbrock Allocation Tests" begin
    # Vector problem - use in-place form for proper cache testing
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
                    vector_prob, solver, dt = 0.1,
                    save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6
                )
                step!(integrator)  # Setup step - initializes caches, may allocate

                # Test perform_step! directly - this is the core stepping function
                # that should be allocation-free (unlike step! which includes saving)
                cache = integrator.cache
                allocs = check_allocs(OrdinaryDiffEqCore.perform_step!,
                    (typeof(integrator), typeof(cache)))

                # These solvers should be allocation-free in perform_step!
                @test_broken length(allocs) == 0

                # We also don't want any dynamic dispatch from this module
                # (currently broken - there are dynamic dispatches to fix)
                @test_broken count(
                    t -> t isa AllocCheck.DynamicDispatch &&
                        any(s -> contains(string(s.file), "Rosenbrock"), t.backtrace),
                    allocs
                ) == 0

                if length(allocs) > 0
                    println("AllocCheck found $(length(allocs)) allocation sites in $(typeof(solver)) perform_step!:")
                    for (i, alloc) in enumerate(allocs[1:min(3, end)])  # Show first 3
                        println("  $i. $alloc")
                    end
                else
                    println("✓ $(typeof(solver)) perform_step! appears allocation-free with AllocCheck")
                end
            end
        end
    end
end
