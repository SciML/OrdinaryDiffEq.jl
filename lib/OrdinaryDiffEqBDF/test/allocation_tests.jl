using Pkg
Pkg.add("AllocCheck")

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
    # Test problem for standard ODE BDF methods
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))

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

    # DAE problem for DAE solvers
    function dae_f!(resid, du, u, p, t)
        resid[1] = -0.5 * u[1] + u[2] - du[1]
        resid[2] = u[1] - u[2] - du[2]
    end
    du0 = zeros(2)
    differential_vars = [true, false]
    dae_prob = DAEProblem(dae_f!, du0, [1.0, 1.0], (0.0, 1.0), differential_vars = differential_vars)

    # Test all exported BDF solvers for allocation-free behavior
    # Standard ODE BDF methods
    bdf_solvers = [ABDF2(), QNDF1(), QBDF1(), QNDF2(), QBDF2(), QNDF(), QBDF(), FBDF(), MEBDF2()]

    # IMEX/Split methods need SplitODEProblem
    imex_solvers = [SBDF(order = 2), SBDF2(), SBDF3(), SBDF4(), IMEXEuler(), IMEXEulerARK()]

    # DAE methods need DAEProblem
    dae_solvers = [DABDF2(), DImplicitEuler(), DFBDF()]

    @testset "BDF Solver Allocation Analysis" begin
        for solver in bdf_solvers
            @testset "$(typeof(solver)) allocation check" begin
                integrator = init(prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                step!(integrator)  # Setup step may allocate

                # Use AllocCheck to verify step! is allocation-free
                allocs = check_allocs(step!, (typeof(integrator),))

                # These solvers should be allocation-free, but mark as broken for now
                # to verify with AllocCheck (more accurate than @allocated)
                @test length(allocs) == 0 broken = true

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
                integrator = init(split_prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                step!(integrator)  # Setup step may allocate

                # Use AllocCheck to verify step! is allocation-free
                allocs = check_allocs(step!, (typeof(integrator),))

                # These solvers should be allocation-free, but mark as broken for now
                # to verify with AllocCheck (more accurate than @allocated)
                @test length(allocs) == 0 broken = true

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

    @testset "DAE Solver Allocation Analysis" begin
        for solver in dae_solvers
            @testset "$(typeof(solver)) allocation check" begin
                integrator = init(dae_prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                step!(integrator)  # Setup step may allocate

                # Use AllocCheck to verify step! is allocation-free
                allocs = check_allocs(step!, (typeof(integrator),))

                # These solvers should be allocation-free, but mark as broken for now
                # to verify with AllocCheck (more accurate than @allocated)
                @test length(allocs) == 0 broken = true

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
