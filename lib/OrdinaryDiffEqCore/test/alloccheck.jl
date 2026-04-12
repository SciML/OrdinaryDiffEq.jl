using OrdinaryDiffEqCore
using OrdinaryDiffEqTsit5
using SciMLBase: FullSpecialize
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqCore infrastructure using AllocCheck.jl.
Tests the core stepping infrastructure functions (loopheader!, _loopfooter!,
apply_step!) independently of any specific algorithm's perform_step!.
These functions are called once per step for every solver, so they must be
allocation-free.

Uses Tsit5 to construct a representative ODEIntegrator.
initdt is excluded per https://github.com/SciML/OrdinaryDiffEq.jl/issues/3231
(intentional allocation that requires cache array refactoring).
"""

@testset "OrdinaryDiffEqCore Infrastructure Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end

    # FullSpecialize avoids FunctionWrappers dynamic dispatch noise
    prob = ODEProblem{true, FullSpecialize}(simple_system!, [1.0, 1.0], (0.0, 10.0))

    integrator = init(
        prob, Tsit5(), dt = 0.1, save_everystep = false,
        abstol = 1.0e-6, reltol = 1.0e-6
    )
    # Warm up: advance so all code paths (accept/reject, controller) are reached
    for _ in 1:5
        step!(integrator)
    end

    # Runtime allocation tests: run FIRST before AllocCheck static analysis,
    # since check_allocs can invalidate compiled code and cause false positives.
    @testset "step! Runtime Allocation Check (save_everystep=false)" begin
        long_integrator = init(
            prob, Tsit5(), dt = 0.1, save_everystep = false,
            abstol = 1.0e-6, reltol = 1.0e-6
        )
        # Warm up to ensure all paths compiled and caches initialized
        for _ in 1:20
            step!(long_integrator)
        end
        GC.gc()
        allocs = 0
        for _ in 1:20
            allocs += @allocated step!(long_integrator)
        end
        # loopheader!/loopfooter!/apply_step! allocate; fix tracked in issue #3231
        @test allocs == 0 broken = true
    end

    # Static allocation analysis for core infrastructure functions.
    # These are called once per step for every ODE solver.

    @testset "loopheader! Static Analysis" begin
        allocs = check_allocs(
            OrdinaryDiffEqCore.loopheader!,
            (typeof(integrator),)
        )
        @test length(allocs) == 0 broken = true
        if length(allocs) > 0
            println("AllocCheck found $(length(allocs)) allocation sites in loopheader!")
        else
            println("loopheader! appears allocation-free with AllocCheck")
        end
    end

    @testset "_loopfooter! Static Analysis" begin
        allocs = check_allocs(
            OrdinaryDiffEqCore._loopfooter!,
            (typeof(integrator),)
        )
        @test length(allocs) == 0 broken = true
        if length(allocs) > 0
            println("AllocCheck found $(length(allocs)) allocation sites in _loopfooter!")
        else
            println("_loopfooter! appears allocation-free with AllocCheck")
        end
    end

    @testset "apply_step! Static Analysis" begin
        allocs = check_allocs(
            OrdinaryDiffEqCore.apply_step!,
            (typeof(integrator),)
        )
        @test length(allocs) == 0 broken = true
        if length(allocs) > 0
            println("AllocCheck found $(length(allocs)) allocation sites in apply_step!")
        else
            println("apply_step! appears allocation-free with AllocCheck")
        end
    end

    @testset "fix_dt_at_bounds! Static Analysis" begin
        allocs = check_allocs(
            OrdinaryDiffEqCore.fix_dt_at_bounds!,
            (typeof(integrator),)
        )
        @test length(allocs) == 0
        if length(allocs) > 0
            println("AllocCheck found $(length(allocs)) allocation sites in fix_dt_at_bounds!")
        else
            println("fix_dt_at_bounds! appears allocation-free with AllocCheck")
        end
    end

    @testset "modify_dt_for_tstops! Static Analysis" begin
        allocs = check_allocs(
            OrdinaryDiffEqCore.modify_dt_for_tstops!,
            (typeof(integrator),)
        )
        @test length(allocs) == 0
        if length(allocs) > 0
            println(
                "AllocCheck found $(length(allocs)) allocation sites in modify_dt_for_tstops!"
            )
        else
            println("modify_dt_for_tstops! appears allocation-free with AllocCheck")
        end
    end

    @testset "perform_step! (Tsit5) Static Analysis" begin
        cache = integrator.cache
        allocs = check_allocs(
            OrdinaryDiffEqCore.perform_step!,
            (typeof(integrator), typeof(cache))
        )
        @test length(allocs) == 0
        if length(allocs) > 0
            println(
                "AllocCheck found $(length(allocs)) allocation sites in Tsit5 perform_step!"
            )
        else
            println("Tsit5 perform_step! appears allocation-free with AllocCheck")
        end
    end
end
