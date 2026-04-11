using Pkg
Pkg.add("AllocCheck")

using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqSDIRK solvers using AllocCheck.jl.
Tests perform_step! directly (the core stepping function) rather than step!,
since step! includes saving operations that naturally allocate.
Uses FullSpecialize to avoid FunctionWrappers dynamic dispatch noise.

All SDIRK solvers are marked broken=true because implicit methods require
linear solves in perform_step!, which currently allocate. Fixing these would
require exposing additional cache arrays for linear solve workspaces.
"""

@testset "SDIRK Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise
    prob = ODEProblem{true, FullSpecialize}(simple_system!, [1.0, 1.0], (0.0, 1.0))

    sdirk_solvers = [
        ImplicitEuler(), ImplicitMidpoint(), Trapezoid(), TRBDF2(), SDIRK2(), SDIRK22(),
        Kvaerno3(), KenCarp3(), Cash4(), Hairer4(), Hairer42(), SSPSDIRK2(),
        Kvaerno4(), Kvaerno5(), KenCarp4(), KenCarp47(), KenCarp5(), KenCarp58(),
        ESDIRK54I8L2SA(), SFSDIRK4(), SFSDIRK5(), CFNLIRK3(),
        SFSDIRK6(), SFSDIRK7(), SFSDIRK8(),
        ESDIRK436L2SA2(), ESDIRK437L2SA(), ESDIRK547L2SA2(), ESDIRK659L2SA(),
    ]

    @testset "SDIRK perform_step! Static Analysis" begin
        for solver in sdirk_solvers
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
