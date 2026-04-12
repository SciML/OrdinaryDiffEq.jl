using OrdinaryDiffEqMultirate
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize, SplitFunction, ODEFunction
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqMultirate solvers using AllocCheck.jl.
Tests perform_step! directly (the core stepping function) rather than step!,
since step! includes saving operations that naturally allocate.
Uses FullSpecialize to avoid FunctionWrappers dynamic dispatch noise.

MREEF is a multirate extrapolated explicit Euler method that uses a
SplitODEProblem (fast + slow components). It is marked broken=true as its
multirate substep structure currently allocates in perform_step!.
"""

@testset "Multirate Allocation Tests" begin
    function f_fast!(du, u, p, t)
        du[1] = -0.9 * u[1]
        du[2] = -0.9 * u[2]
    end
    function f_slow!(du, u, p, t)
        du[1] = -0.1 * u[1]
        du[2] = -0.1 * u[2]
    end

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise
    prob = SplitODEProblem(
        SplitFunction(
            ODEFunction{true, FullSpecialize}(f_fast!),
            ODEFunction{true, FullSpecialize}(f_slow!)
        ),
        [1.0, 1.0], (0.0, 1.0)
    )

    mreef_solvers = [MREEF(), MREEF(m = 8, order = 3)]

    @testset "MREEF perform_step! Static Analysis" begin
        for solver in mreef_solvers
            @testset "$(solver) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.1, save_everystep = false, adaptive = false
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
                        "AllocCheck found $(length(allocs)) allocation sites in MREEF perform_step!"
                    )
                else
                    println("MREEF perform_step! appears allocation-free with AllocCheck")
                end
            end
        end
    end
end
