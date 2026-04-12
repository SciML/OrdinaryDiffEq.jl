using OrdinaryDiffEqNewmark
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize, DynamicalODEFunction
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqNewmark solvers using AllocCheck.jl.
Tests perform_step! directly (the core stepping function) rather than step!,
since step! includes saving operations that naturally allocate.

NewmarkBeta is an implicit structural dynamics method marked broken=true
because it requires a linear solve in perform_step! that currently allocates.
"""

@testset "Newmark Allocation Tests" begin
    # Harmonic oscillator: dv/dt = -u, du/dt = v
    function f1!(dv, v, u, p, t)
        dv .= -u
    end
    function f2!(du, v, u, p, t)
        du .= v
    end

    v0 = [1.0, 0.0]
    u0 = [0.0, 1.0]

    # FullSpecialize on DynamicalODEFunction avoids FunctionWrappers noise
    ff = DynamicalODEFunction{true, FullSpecialize}(f1!, f2!)
    prob = DynamicalODEProblem(ff, v0, u0, (0.0, 1.0))

    @testset "NewmarkBeta perform_step! Static Analysis" begin
        integrator = init(
            prob, NewmarkBeta(), dt = 0.1, save_everystep = false, adaptive = false
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
                "AllocCheck found $(length(allocs)) allocation sites in NewmarkBeta perform_step!"
            )
        else
            println("NewmarkBeta perform_step! appears allocation-free with AllocCheck")
        end
    end
end
