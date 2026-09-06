using OrdinaryDiffEqPDIRK, OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NLNewton, NLFunctional, NLAnderson,
    NonlinearSolveAlg, HomotopyNonlinearSolveAlg
using NonlinearSolve: NewtonRaphson
using Test
import OrdinaryDiffEqCore

@test Threads.nthreads() >= 4

stats_counters(stats) = NamedTuple{fieldnames(typeof(stats))}(
    Tuple(getfield(stats, field) for field in fieldnames(typeof(stats)))
)

function stats_rhs!(du, u, p, t)
    @. du = -u * u
    return nothing
end
stats_rhs(u, p, t) = -u .* u

@testset "PDIRK concurrency" begin
    @testset "PDIRK statistics agree across threading modes" begin
        for rhs in (stats_rhs!, stats_rhs), nl in (
                    NLNewton(), NLFunctional(), NLAnderson(), NonlinearSolveAlg(NewtonRaphson()),
                    HomotopyNonlinearSolveAlg(),
                )
            @testset "$(nameof(typeof(nl))), $(rhs)" begin
                prob = ODEProblem(rhs, ones(8), (0.0, 1.0))
                serial = solve(
                    prob, PDIRK44(; threading = false, nlsolve = nl);
                    dt = 0.01, adaptive = false, save_everystep = false
                )
                @test successful_retcode(serial)
                expected = stats_counters(serial.stats)
                for _ in 1:30
                    threaded = solve(
                        prob, PDIRK44(; threading = true, nlsolve = nl);
                        dt = 0.01, adaptive = false, save_everystep = false
                    )
                    @test successful_retcode(threaded)
                    @test threaded.u[end] ≈ serial.u[end]
                    @test stats_counters(threaded.stats) == expected
                end
            end
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
        function rhs(u, p, t)
            yield()
            return -u .* u
        end
        function rhs!(du, u, p, t)
            yield()
            @. du = -u * u
            return nothing
        end
        for f in (rhs!, rhs)
            prob = ODEProblem(f, ones(8), (0.0, 1.0))
            serial = solve(
                prob, PDIRK44(; threading = false, nlsolve = NLNewton(; always_new = true));
                dt = 0.01, adaptive = false
            )
            @test successful_retcode(serial)
            expected = stats_counters(serial.stats)
            for _ in 1:30
                threaded = solve(
                    prob, PDIRK44(; threading = true, nlsolve = NLNewton(; always_new = true));
                    dt = 0.01, adaptive = false
                )
                @test stats_counters(threaded.stats) == expected
                @test successful_retcode(threaded)
            end
        end
    end
end
