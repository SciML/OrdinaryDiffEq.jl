using OrdinaryDiffEqExtrapolation
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqExtrapolation solvers using AllocCheck.jl.
Tests perform_step! directly (the core stepping function) rather than step!,
since step! includes saving operations that naturally allocate.
Uses FullSpecialize to avoid FunctionWrappers dynamic dispatch noise.

All extrapolation methods are marked broken=true. These methods adaptively
select the extrapolation order and maintain multiple substep solutions,
requiring dynamic allocation in perform_step!. Fixing would require
pre-allocating all substep work arrays in the cache.
"""

@testset "Extrapolation Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise
    prob = ODEProblem{true, FullSpecialize}(simple_system!, [1.0, 1.0], (0.0, 1.0))

    explicit_extrap_solvers = [
        AitkenNeville(), ExtrapolationMidpointDeuflhard(),
        ExtrapolationMidpointHairerWanner(),
    ]

    implicit_extrap_solvers = [
        ImplicitEulerExtrapolation(), ImplicitDeuflhardExtrapolation(),
        ImplicitHairerWannerExtrapolation(), ImplicitEulerBarycentricExtrapolation(),
    ]

    @testset "Extrapolation perform_step! Static Analysis" begin
        for solver in explicit_extrap_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.1, save_everystep = false,
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
                        "AllocCheck found $(length(allocs)) allocation sites in $(typeof(solver)) perform_step!"
                    )
                else
                    println(
                        "$(typeof(solver)) perform_step! appears allocation-free with AllocCheck"
                    )
                end
            end
        end

        for solver in implicit_extrap_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.1, save_everystep = false,
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
