using OrdinaryDiffEqFunctionMap
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize, DiscreteFunction
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqFunctionMap (FunctionMap) using AllocCheck.jl.
Tests perform_step! directly (the core stepping function) rather than step!,
since step! includes saving operations that naturally allocate.
FunctionMap solves DiscreteProblem: u[n+1] = f(u[n], p, t).
Uses DiscreteFunction{true, FullSpecialize} to avoid FunctionWrappers noise.
(DiscreteProblem does not accept {true, FullSpecialize} type parameters directly.)
"""

@testset "FunctionMap Allocation Tests" begin
    function discrete_map!(du, u, p, t)
        du[1] = 0.9 * u[1]
        du[2] = 0.8 * u[2]
    end

    # DiscreteProblem does not support {true, FullSpecialize} type parameters;
    # use DiscreteFunction{true, FullSpecialize} wrapper instead.
    ff = DiscreteFunction{true, FullSpecialize}(discrete_map!)
    prob = DiscreteProblem(ff, [1.0, 1.0], (0.0, 10.0))

    @testset "FunctionMap perform_step! Static Analysis" begin
        integrator = init(prob, FunctionMap(), save_everystep = false)
        step!(integrator)

        cache = integrator.cache
        allocs = check_allocs(
            OrdinaryDiffEqCore.perform_step!,
            (typeof(integrator), typeof(cache))
        )

        @test length(allocs) == 0

        if length(allocs) > 0
            println(
                "AllocCheck found $(length(allocs)) allocation sites in FunctionMap perform_step!"
            )
        else
            println("FunctionMap perform_step! appears allocation-free with AllocCheck")
        end
    end

    @testset "FunctionMap(scale_by_time=true) perform_step! Static Analysis" begin
        integrator = init(
            prob, FunctionMap(scale_by_time = true), save_everystep = false
        )
        step!(integrator)

        cache = integrator.cache
        allocs = check_allocs(
            OrdinaryDiffEqCore.perform_step!,
            (typeof(integrator), typeof(cache))
        )

        @test length(allocs) == 0

        if length(allocs) > 0
            println(
                "AllocCheck found $(length(allocs)) allocation sites in FunctionMap(scale_by_time=true) perform_step!"
            )
        else
            println(
                "FunctionMap(scale_by_time=true) perform_step! appears allocation-free with AllocCheck"
            )
        end
    end
end
