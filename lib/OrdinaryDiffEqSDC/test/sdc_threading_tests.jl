using OrdinaryDiffEqSDC
using SciMLBase
using Test

# A diagonal QΔ removes the only node-to-node coupling inside a sweep, so the M
# node solves become independent and can run concurrently. The results must not
# depend on whether they did.
rotation!(du, u, p, t) = (du[1] = -u[2]; du[2] = u[1]; nothing)
const ROTATION = ODEProblem(rotation!, [1.0, 0.0], (0.0, 2 * π))

sdc(sweeper, threading; M = 4, K = 6) = SDC(
    num_nodes = M, num_sweeps = K, sweeper = sweeper, threading = threading
)

@testset "threading is rejected for sweepers that couple the nodes" begin
    for sweeper in (SDCSweeper.BE, SDCSweeper.FE, SDCSweeper.Trapezoid, SDCSweeper.LU)
        @test_throws ArgumentError sdc(sweeper, true)
        @test sdc(sweeper, false) isa SDC
    end
    for sweeper in OrdinaryDiffEqSDC.SDC_DIAGONAL_SWEEPERS
        @test sdc(sweeper, true) isa SDC
    end
end

@testset "threaded sweeps reproduce the sequential result exactly" begin
    # Bitwise, not approximately: each node writes only its own slot, so
    # concurrency cannot change the arithmetic or its order.
    for sweeper in OrdinaryDiffEqSDC.SDC_DIAGONAL_SWEEPERS
        sequential = solve(
            ROTATION, sdc(sweeper, false); dt = (2 * π) / 64, adaptive = false
        )
        threaded = solve(
            ROTATION, sdc(sweeper, true); dt = (2 * π) / 64, adaptive = false
        )
        # Statistics are deliberately not compared. `nlsolve!` and `calc_W!`
        # accumulate `integrator.stats` with a plain `+=` from whichever thread
        # runs the node, so the counters undercount when threaded.
        @testset "$sweeper" begin
            @test threaded.u[end] == sequential.u[end]
        end
    end
end

@testset "threaded sweeps work with adaptivity and out-of-place problems" begin
    alg = sdc(SDCSweeper.MIN_SR_S, true)
    adaptive = solve(ROTATION, alg; abstol = 1.0e-9, reltol = 1.0e-9)
    @test SciMLBase.successful_retcode(adaptive)
    @test maximum(abs.(adaptive.u[end] .- [cos(2 * π), sin(2 * π)])) < 1.0e-7

    oop = ODEProblem((u, p, t) -> [-u[2], u[1]], [1.0, 0.0], (0.0, 2 * π))
    sequential = solve(oop, sdc(SDCSweeper.MIN_SR_S, false); dt = 0.1, adaptive = false)
    threaded = solve(oop, sdc(SDCSweeper.MIN_SR_S, true); dt = 0.1, adaptive = false)
    @test threaded.u[end] == sequential.u[end]
end

@testset "a failed node solve is not lost between threads" begin
    # `nlsolve!` assigns `integrator.force_stepfail` from every node, so a node
    # that converges after one that did not would otherwise erase the failure.
    # The per-node flags are what the step actually acts on.
    cache_type = typeof(
        init(
            ROTATION, sdc(SDCSweeper.MIN_SR_S, true);
            dt = 0.1, adaptive = false
        ).cache
    )
    @test hasfield(cache_type, :failed)

    integrator = init(
        ROTATION, sdc(SDCSweeper.MIN_SR_S, true); dt = 0.1, adaptive = false
    )
    step!(integrator)
    @test !any(integrator.cache.failed)
    @test !integrator.force_stepfail
end
