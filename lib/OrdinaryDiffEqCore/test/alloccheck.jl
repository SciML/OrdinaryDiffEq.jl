using OrdinaryDiffEqCore
using OrdinaryDiffEqTsit5
using SciMLBase: FullSpecialize
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqCore infrastructure using AllocCheck.jl.

Tests the core stepping infrastructure functions called once per step for every solver:
  - loopheader!
  - _loopfooter!
  - apply_step!
  - fix_dt_at_bounds!
  - modify_dt_for_tstops!
  - perform_step! (via Tsit5 as the reference solver)

initdt is excluded: it currently allocates intentionally and is tracked separately
in https://github.com/SciML/OrdinaryDiffEq.jl/issues/3231.

Uses FullSpecialize to avoid FunctionWrappers false-positive allocations.
Runtime @allocated tests run FIRST before check_allocs static analysis,
since check_allocs can invalidate compiled code and cause spurious allocations.
"""

@testset "OrdinaryDiffEqCore Infrastructure AllocCheck" begin
    function f!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end

    prob = ODEProblem{true, FullSpecialize}(f!, [1.0, 1.0], (0.0, 10.0))
    integrator = init(
        prob, Tsit5(), dt = 0.1, save_everystep = false,
        abstol = 1e-6, reltol = 1e-6
    )
    # Warm up: advance enough steps so all code paths (accept, controller) are exercised
    for _ in 1:10
        step!(integrator)
    end

    # --- Runtime allocation checks (must run before check_allocs) ---

    @testset "step! runtime allocation check (save_everystep=false)" begin
        long_integrator = init(
            prob, Tsit5(), dt = 0.1, save_everystep = false,
            abstol = 1e-6, reltol = 1e-6
        )
        for _ in 1:20
            step!(long_integrator)
        end
        GC.gc()
        allocs = sum(@allocated(step!(long_integrator)) for _ in 1:20)
        @test allocs == 0 broken = true
    end

    # --- Static allocation analysis (check_allocs) ---

    @testset "loopheader! alloccheck" begin
        allocs = check_allocs(OrdinaryDiffEqCore.loopheader!, (typeof(integrator),))
        @test length(allocs) == 0 broken = true
        if !isempty(allocs)
            println("loopheader! has $(length(allocs)) allocation site(s)")
        end
    end

    @testset "_loopfooter! alloccheck" begin
        allocs = check_allocs(OrdinaryDiffEqCore._loopfooter!, (typeof(integrator),))
        @test length(allocs) == 0 broken = true
        if !isempty(allocs)
            println("_loopfooter! has $(length(allocs)) allocation site(s)")
        end
    end

    @testset "apply_step! alloccheck" begin
        allocs = check_allocs(OrdinaryDiffEqCore.apply_step!, (typeof(integrator),))
        @test length(allocs) == 0 broken = true
        if !isempty(allocs)
            println("apply_step! has $(length(allocs)) allocation site(s)")
        end
    end

    @testset "fix_dt_at_bounds! alloccheck" begin
        allocs = check_allocs(OrdinaryDiffEqCore.fix_dt_at_bounds!, (typeof(integrator),))
        @test length(allocs) == 0 broken = true
        if !isempty(allocs)
            println("fix_dt_at_bounds! has $(length(allocs)) allocation site(s)")
        end
    end

    @testset "modify_dt_for_tstops! alloccheck" begin
        allocs = check_allocs(
            OrdinaryDiffEqCore.modify_dt_for_tstops!, (typeof(integrator),))
        @test length(allocs) == 0 broken = true
        if !isempty(allocs)
            println("modify_dt_for_tstops! has $(length(allocs)) allocation site(s)")
        end
    end

    @testset "perform_step! (Tsit5) alloccheck" begin
        cache = integrator.cache
        allocs = check_allocs(
            OrdinaryDiffEqCore.perform_step!,
            (typeof(integrator), typeof(cache))
        )
        @test length(allocs) == 0
        if !isempty(allocs)
            println("Tsit5 perform_step! has $(length(allocs)) allocation site(s)")
        end
    end
end
