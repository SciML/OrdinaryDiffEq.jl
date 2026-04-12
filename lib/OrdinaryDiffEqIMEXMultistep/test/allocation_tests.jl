using OrdinaryDiffEqIMEXMultistep
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize, SplitFunction, ODEFunction
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqIMEXMultistep solvers using AllocCheck.jl.
Tests perform_step! directly (the core stepping function) rather than step!,
since step! includes saving operations that naturally allocate.
Uses FullSpecialize to avoid FunctionWrappers dynamic dispatch noise.

CNAB2 and CNLF2 are IMEX multistep methods and are marked broken=true
because their implicit part requires linear solves that currently allocate.
"""

@testset "IMEXMultistep Allocation Tests" begin
    function f1!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = 0.0
    end
    function f2!(du, u, p, t)
        du[1] = 0.0
        du[2] = -1.5 * u[2]
    end

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise
    prob = SplitODEProblem(
        SplitFunction(
            ODEFunction{true, FullSpecialize}(f1!),
            ODEFunction{true, FullSpecialize}(f2!)
        ),
        [1.0, 1.0], (0.0, 1.0)
    )

    imex_ms_solvers = [CNAB2(), CNLF2()]

    @testset "IMEXMultistep perform_step! Static Analysis" begin
        for solver in imex_ms_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.1, save_everystep = false,
                    abstol = 1.0e-6, reltol = 1.0e-6
                )
                # Multistep: advance past startup phase
                for _ in 1:3
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
