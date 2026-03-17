using Test
using StochasticDiffEq
using Random
using AllocCheck

"""
    Allocation tests for StochasticDiffEq.jl

These tests ensure that the hot path functions (perform_step!, loopheader!, loopfooter!)
don't allocate memory during stepping, which is critical for performance.
"""

# Wrapper that discards the return value of step!.
# ODE's _step! returns integrator.sol.retcode, and on Julia 1.10 boxing that
# ReturnCode.T return value inside @allocated costs 16 bytes. Discarding it
# avoids measuring return-value boxing as a false positive.
@inline function step_void!(integrator)
    step!(integrator)
    return nothing
end

@testset "Allocation Tests" begin
    Random.seed!(12345)

    # Test problem definitions (in-place, 3D diagonal noise)
    f_iip!(du, u, p, t) = (du[1] = 0.1 * u[1]; du[2] = 0.2 * u[2]; du[3] = 0.3 * u[3])
    g_iip!(du, u, p, t) = (du[1] = 0.2 * u[1]; du[2] = 0.3 * u[2]; du[3] = 0.4 * u[3])
    u0 = [1.0, 1.0, 1.0]
    tspan = (0.0, 10.0)  # Long enough to ensure many steps
    prob_iip = SDEProblem(f_iip!, g_iip!, u0, tspan)

    @testset "EM stepping allocation" begin
        # Initialize integrator with save_on=false to avoid solution saving allocations
        integrator = init(prob_iip, EM(), dt = 0.01, save_on = false)

        # Warm up to ensure JIT compilation
        for _ in 1:10
            step_void!(integrator)
        end

        # Test allocation per step
        allocs_per_step = @allocated step_void!(integrator)
        @test allocs_per_step == 0
    end

    @testset "EulerHeun stepping allocation" begin
        integrator = init(prob_iip, EulerHeun(), dt = 0.01, save_on = false)

        # Warm up
        for _ in 1:10
            step_void!(integrator)
        end

        allocs_per_step = @allocated step_void!(integrator)
        @test allocs_per_step == 0
    end

    @testset "RKMil stepping allocation" begin
        integrator = init(prob_iip, RKMil(), dt = 0.01, save_on = false)

        # Warm up
        for _ in 1:10
            step_void!(integrator)
        end

        allocs_per_step = @allocated step_void!(integrator)
        @test allocs_per_step == 0
    end

    @testset "SOSRI stepping allocation" begin
        integrator = init(prob_iip, SOSRI(), save_on = false)

        # Warm up
        for _ in 1:10
            step_void!(integrator)
        end

        allocs_per_step = @allocated step_void!(integrator)
        @test allocs_per_step == 0
    end

    @testset "SOSRA stepping allocation" begin
        integrator = init(prob_iip, SOSRA(), save_on = false)

        # Warm up
        for _ in 1:10
            step_void!(integrator)
        end

        allocs_per_step = @allocated step_void!(integrator)
        @test allocs_per_step == 0
    end

    @testset "SKenCarp stepping allocation" begin
        integrator = init(prob_iip, SKenCarp(), dt = 0.01, adaptive = false, save_on = false)

        # SKenCarp has a deep call graph (NL solver, Jacobian, linear solve)
        # that needs extra warmup to fully populate method dispatch caches
        for _ in 1:50
            step_void!(integrator)
        end

        allocs_per_step = minimum(@allocated(step_void!(integrator)) for _ in 1:5)
        @test allocs_per_step == 0
    end

    # Test with scalar SDE (out-of-place)
    @testset "Scalar SDE allocations" begin
        f_scalar(u, p, t) = 0.1 * u
        g_scalar(u, p, t) = 0.2 * u
        prob_scalar = SDEProblem(f_scalar, g_scalar, 1.0, tspan)

        integrator = init(prob_scalar, EM(), dt = 0.01, save_on = false)

        # Warm up
        for _ in 1:10
            step_void!(integrator)
        end

        allocs_per_step = @allocated step_void!(integrator)
        @test allocs_per_step == 0
    end

    # Test perform_step! directly
    @testset "perform_step! allocation - EM" begin
        integrator = init(prob_iip, EM(), dt = 0.01, save_on = false)
        cache = integrator.cache

        # Warm up
        for _ in 1:10
            StochasticDiffEq.perform_step!(integrator, cache)
            # Manually update t to avoid going past tspan
            integrator.t = mod(integrator.t, 9.0)
        end

        allocs = @allocated StochasticDiffEq.perform_step!(integrator, cache)
        @test allocs == 0
    end

    @testset "loopheader! allocation" begin
        integrator = init(prob_iip, EM(), dt = 0.01, save_on = false)

        # Warm up
        for _ in 1:10
            StochasticDiffEq.loopheader!(integrator)
        end

        allocs = @allocated StochasticDiffEq.loopheader!(integrator)
        @test allocs == 0
    end

    @testset "loopfooter! allocation" begin
        integrator = init(prob_iip, EM(), dt = 0.01, save_on = false)

        # Warm up
        for _ in 1:10
            StochasticDiffEq.loopfooter!(integrator)
        end

        allocs = @allocated StochasticDiffEq.loopfooter!(integrator)
        @test allocs == 0
    end
end
