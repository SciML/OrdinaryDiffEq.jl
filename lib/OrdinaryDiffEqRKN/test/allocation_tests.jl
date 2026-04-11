using Pkg
Pkg.add("AllocCheck")

using OrdinaryDiffEqRKN
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize, DynamicalODEFunction
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqRKN solvers using AllocCheck.jl.
Tests perform_step! directly (the core stepping function) rather than step!,
since step! includes saving operations that naturally allocate.
Uses FullSpecialize on DynamicalODEFunction to avoid FunctionWrappers noise.

RKN methods solve second-order ODEs in the form:
  dv/dt = f1(v, u, p, t)   (acceleration/force)
  du/dt = f2(v, u, p, t)   (velocity, often f2 = v)

All RKN methods are marked broken=true. The DynamicalODEProblem perform_step!
involves ArrayPartition construction per step that currently allocates.
Fixing requires refactoring the partitioned-state update in perform_step!.
"""

@testset "RKN Allocation Tests" begin
    # Simple harmonic oscillator: dv = -u, du = v
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

    # All RKN methods — both explicit and implicit currently allocate
    all_rkn_solvers = [
        Nystrom4(), FineRKN4(), FineRKN5(),
        Nystrom4VelocityIndependent(), Nystrom5VelocityIndependent(),
        DPRKN4(), DPRKN5(), DPRKN6(), DPRKN6FM(), DPRKN8(), DPRKN12(),
        ERKN4(), ERKN5(), ERKN7(), RKN4(),
        IRKN3(), IRKN4(),
    ]

    @testset "RKN perform_step! Static Analysis" begin
        for solver in all_rkn_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.01, save_everystep = false,
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
