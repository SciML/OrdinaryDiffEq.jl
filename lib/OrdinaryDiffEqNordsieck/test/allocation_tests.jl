using OrdinaryDiffEqNordsieck
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqNordsieck solvers using AllocCheck.jl.
Tests perform_step! directly (the core stepping function) rather than step!,
since step! includes saving operations that naturally allocate.
Uses FullSpecialize to avoid FunctionWrappers dynamic dispatch noise.

All Nordsieck/VODE methods are marked broken=true. These variable-order,
variable-step methods manage a Nordsieck array and dynamically resize
internal history structures, which currently allocate in perform_step!.
"""

@testset "Nordsieck Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise
    prob = ODEProblem{true, FullSpecialize}(simple_system!, [1.0, 1.0], (0.0, 10.0))

    nordsieck_solvers = [AN5(), JVODE_Adams(), JVODE_BDF()]

    @testset "Nordsieck perform_step! Static Analysis" begin
        for solver in nordsieck_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.1, save_everystep = false,
                    abstol = 1.0e-6, reltol = 1.0e-6
                )
                for _ in 1:5
                    step!(integrator)
                end

                cache = integrator.cache
                allocs = check_allocs(
                    OrdinaryDiffEqCore.perform_step!,
                    (typeof(integrator), typeof(cache))
                )

                @test length(allocs) == 0 broken = true

                if length(allocs) > 0
                    println(
                        "AllocCheck found $(length(allocs)) allocation sites in $(typeof(solver)) perform_step!"
                    )
                else
                    println(
                        "$(typeof(solver)) perform_step! appears allocation-free with AllocCheck"
                    )
                end
            end
        end
    end
end
