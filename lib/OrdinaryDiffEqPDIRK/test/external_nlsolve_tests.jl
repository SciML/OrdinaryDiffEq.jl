module ExternalNLSolverFactoryClient

    using OrdinaryDiffEqCore: AbstractNLSolverAlgorithm
    using OrdinaryDiffEqNonlinearSolve: NLNewton
    import OrdinaryDiffEqNonlinearSolve: build_nlsolver

    struct PositionalNewton <: AbstractNLSolverAlgorithm end
    struct BufferedNewton <: AbstractNLSolverAlgorithm end

    build_nlsolver(alg, ::PositionalNewton, args...) =
        build_nlsolver(alg, NLNewton(), args...)

    function build_nlsolver(alg, ::BufferedNewton, args...; stats_delta = nothing)
        return build_nlsolver(alg, NLNewton(), args...; stats_delta)
    end

end

using .ExternalNLSolverFactoryClient: PositionalNewton, BufferedNewton
using OrdinaryDiffEqPDIRK
using Test

@testset "External nonlinear-solver factories" begin
    oop = ODEProblem((u, p, t) -> -u, [1.0], (0.0, 1.0))
    iip = ODEProblem((du, u, p, t) -> (du .= -u), [1.0], (0.0, 1.0))
    for (inplace, prob) in ((false, oop), (true, iip))
        reference = solve(
            prob, PDIRK44(threading = false); dt = 0.01, adaptive = false
        )
        @test successful_retcode(reference)
        @testset "Positional factory, $(inplace)" begin
            sol = solve(
                prob, PDIRK44(threading = false, nlsolve = PositionalNewton());
                dt = 0.01, adaptive = false
            )
            @test successful_retcode(sol)
            @test sol.u[end] ≈ [exp(-1)] atol = 1.0e-8
            @test sol.u == reference.u
            for field in fieldnames(typeof(reference.stats))
                @test getfield(sol.stats, field) == getfield(reference.stats, field)
            end
            @test_throws MethodError init(
                prob, PDIRK44(threading = true, nlsolve = PositionalNewton()); dt = 0.01
            )
        end
        @testset "Buffer-aware factory, $(inplace)" begin
            for threading in (false, true)
                sol = solve(
                    prob, PDIRK44(; threading, nlsolve = BufferedNewton());
                    dt = 0.01, adaptive = false
                )
                @test successful_retcode(sol)
                @test sol.u == reference.u
                for field in fieldnames(typeof(reference.stats))
                    @test getfield(sol.stats, field) == getfield(reference.stats, field)
                end
            end
        end
    end
end
