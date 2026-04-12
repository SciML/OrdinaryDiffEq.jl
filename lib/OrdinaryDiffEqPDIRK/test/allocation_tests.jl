using OrdinaryDiffEqPDIRK
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqPDIRK solvers using AllocCheck.jl.
Tests perform_step! directly (the core stepping function) rather than step!,
since step! includes saving operations that naturally allocate.
Uses FullSpecialize to avoid FunctionWrappers dynamic dispatch noise.

PDIRK44 is an implicit parallel diagonally-implicit RK method and is marked
broken=true because implicit methods require linear solves that currently allocate.
"""

@testset "PDIRK Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise
    prob = ODEProblem{true, FullSpecialize}(simple_system!, [1.0, 1.0], (0.0, 1.0))

    @testset "PDIRK44 perform_step! Static Analysis" begin
        integrator = init(
            prob, PDIRK44(), dt = 0.1, save_everystep = false,
            abstol = 1.0e-6, reltol = 1.0e-6
        )
        step!(integrator)

        cache = integrator.cache
        allocs = check_allocs(
            OrdinaryDiffEqCore.perform_step!,
            (typeof(integrator), typeof(cache))
        )

        @test length(allocs) == 0 broken = true

        if length(allocs) > 0
            println(
                "AllocCheck found $(length(allocs)) allocation sites in PDIRK44 perform_step!"
            )
        else
            println("PDIRK44 perform_step! appears allocation-free with AllocCheck")
        end
    end
end
