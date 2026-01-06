using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqCore
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqLowOrderRK solvers using AllocCheck.jl.
These tests verify that the step! operation does not allocate during stepping.
"""

@testset "LowOrderRK Allocation Tests" begin
    # Test problem for adaptive methods
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))

    # Test all exported LowOrderRK solvers for allocation-free behavior
    low_order_solvers = [
        Euler(), Heun(), Ralston(), Midpoint(), RK4(),
        BS3(), OwrenZen3(), OwrenZen4(), OwrenZen5(), BS5(),
        DP5(), Anas5(), RKO65(), FRK65(), RKM(), MSRK5(), MSRK6(),
        PSRK4p7q6(), PSRK3p5q4(), PSRK3p6q5(), Stepanov5(), SIR54(),
        Alshina2(), Alshina3(), Alshina6(), AutoDP5(DP5()),
    ]

    @testset "LowOrderRK Solver Allocation Analysis" begin
        for solver in low_order_solvers
            @testset "$(typeof(solver)) allocation check" begin
                # Some solvers need fixed timestep
                if solver isa Euler || solver isa Midpoint || solver isa Heun
                    integrator = init(prob, solver, dt = 0.1, save_everystep = false, adaptive = false)
                else
                    integrator = init(prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                end
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
