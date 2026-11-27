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
using OrdinaryDiffEqSDIRK: TRBDF2
using Test

@testset "External nonlinear-solver factories" begin
    oop = ODEProblem((u, p, t) -> -u, [1.0], (0.0, 1.0))
    iip = ODEProblem((du, u, p, t) -> (du .= -u), [1.0], (0.0, 1.0))
    for (inplace, prob) in ((false, oop), (true, iip)),
            alg in (PDIRK44(threading = false), PDIRK44(threading = true), TRBDF2())
        reference = solve(prob, alg; dt = 0.01, adaptive = false)
        @test successful_retcode(reference)
        for factory in (PositionalNewton(), BufferedNewton())
            @testset "$(typeof(alg)), $(typeof(factory)), $(inplace)" begin
                custom_alg = alg isa PDIRK44 ?
                    PDIRK44(threading = alg.threading, nlsolve = factory) :
                    TRBDF2(nlsolve = factory)
                sol = solve(prob, custom_alg; dt = 0.01, adaptive = false)
                @test successful_retcode(sol)
                @test sol.u[end] ≈ [exp(-1)] atol = 2.0e-6
                @test sol.u == reference.u
                for field in fieldnames(typeof(reference.stats))
                    @test getfield(sol.stats, field) == getfield(reference.stats, field)
                end
            end
        end
    end
end
