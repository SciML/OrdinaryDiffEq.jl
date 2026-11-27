using OrdinaryDiffEqPDIRK, OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NLNewton, NLFunctional, NLAnderson,
    NonlinearSolveAlg, HomotopyNonlinearSolveAlg
using NonlinearSolve: NewtonRaphson
using LinearSolve: KrylovJL_GMRES
using Test
import OrdinaryDiffEqCore

@test Threads.nthreads() >= 4

using Serialization
using Polyester
using OrdinaryDiffEqCore: BaseThreads, PolyesterThreads
include("threaded_stats_reference.jl")

# One thread runs the two-cache schedule without concurrent counter updates.
references = mktemp() do filename, io
    close(io)
    reference_file = joinpath(@__DIR__, "threaded_stats_reference.jl")
    script = "using Serialization; include(ARGS[1]); serialize(ARGS[2], threaded_stats_reference())"
    run(`$(Base.julia_cmd()) --threads=1 --project=$(dirname(Base.active_project())) -e $script $reference_file $filename`)
    return deserialize(filename)
end

@testset "PDIRK concurrency" begin
    @testset "PDIRK statistics agree across threading modes" begin
        for rhs in (stats_rhs!, stats_rhs), nl in stats_nonlinear_solvers
            @testset "$(nameof(typeof(nl))), $(rhs)" begin
                prob = ODEProblem(rhs, ones(8), (0.0, 1.0))
                serial = solve(
                    prob, PDIRK44(; threading = false, nlsolve = nl);
                    dt = 0.01, adaptive = false, save_everystep = false
                )
                @test successful_retcode(serial)
                expected = references[(nameof(rhs), nameof(typeof(nl)))]
                for mode in (true, BaseThreads(), PolyesterThreads()), _ in 1:30
                    threaded = solve(
                        prob, PDIRK44(; threading = mode, nlsolve = nl);
                        dt = 0.01, adaptive = false, save_everystep = false
                    )
                    @test successful_retcode(threaded)
                    @test threaded.u[end] ≈ serial.u[end]
                    @test stats_counters(threaded.stats) == expected
                end
            end
        end
    end

    @testset "Jacobian-free Krylov statistics agree across threading modes" begin
        prob = ODEProblem(stats_rhs!, ones(8), (0.0, 1.0))
        krylov(threading) = PDIRK44(;
            threading, linsolve = KrylovJL_GMRES(), concrete_jac = false
        )
        serial = solve(
            prob, krylov(false); dt = 0.01, adaptive = false, save_everystep = false
        )
        @test successful_retcode(serial)
        expected = references[:krylov]
        for mode in (true, BaseThreads(), PolyesterThreads()), _ in 1:30
            threaded = solve(
                prob, krylov(mode); dt = 0.01, adaptive = false, save_everystep = false
            )
            @test successful_retcode(threaded)
            @test threaded.u[end] ≈ serial.u[end]
            @test stats_counters(threaded.stats) == expected
        end
    end

    @testset "failed workers merge before PDIRK returns" begin
        for inplace in (false, true), failure_time in (0.005, 0.01 * (2 / 3), 0.01 * (1 / 3))
            rhs(u, p, t) = fill(t == p ? NaN : 0.0, size(u))
            function rhs!(du, u, p, t)
                fill!(du, t == p ? NaN : 0.0)
                return nothing
            end
            prob = ODEProblem(inplace ? rhs! : rhs, ones(8), (0.0, 1.0), failure_time)
            for _ in 1:30
                integrator = init(
                    prob, PDIRK44(; threading = true, nlsolve = NLFunctional());
                    dt = 0.01, adaptive = false
                )
                OrdinaryDiffEqCore.perform_step!(integrator, integrator.cache)
                @test integrator.force_stepfail
                @test count(OrdinaryDiffEqNonlinearSolve.nlsolvefail, integrator.cache.nlsolver) == 1
                @test integrator.stats.nnonlinconvfail == 1
                @test integrator.stats.nf == integrator.stats.nnonliniter > 0
                for solver in integrator.cache.nlsolver
                    sink = OrdinaryDiffEqCore.stats_sink(solver.cache)
                    @test all(iszero(getfield(sink, field)) for field in fieldnames(typeof(sink)))
                end
            end
        end
    end


    @testset "always-new Newton keeps worker Jacobian state private" begin
        for f in (yielding_stats_rhs!, yielding_stats_rhs)
            prob = ODEProblem(f, ones(8), (0.0, 1.0))
            serial = solve(
                prob, PDIRK44(; threading = false, nlsolve = NLNewton(; always_new = true));
                dt = 0.01, adaptive = false
            )
            @test successful_retcode(serial)
            expected = references[nameof(f)]
            for _ in 1:30
                threaded = solve(
                    prob, PDIRK44(; threading = true, nlsolve = NLNewton(; always_new = true));
                    dt = 0.01, adaptive = false
                )
                @test stats_counters(threaded.stats) == expected
                @test threaded.u[end] ≈ serial.u[end]
                @test successful_retcode(threaded)
            end
        end
    end
end
