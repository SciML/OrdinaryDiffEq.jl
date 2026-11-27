using OrdinaryDiffEqPDIRK
using OrdinaryDiffEqPDIRK: Sequential
using Test

@testset "Serial PDIRK work schedule" begin
    function rhs!(du, u, calls, t)
        calls[] += 1
        for i in eachindex(u)
            du[i] = -i * u[i] * (1 + 0.1u[i])
        end
        return nothing
    end
    for threading in (false, Sequential())
        calls = Ref(0)
        prob = ODEProblem(rhs!, ones(50), (0.0, 1.0), calls)
        integrator = init(prob, PDIRK44(; threading); dt = 1 / 200, adaptive = false)
        @test length(integrator.cache.nlsolver) == 1
        sol = solve!(integrator)
        @test successful_retcode(sol)
        @test sol.u[end] ≈ [exp(-i) / (1.1 - 0.1exp(-i)) for i in 1:50] atol = 1.0e-9
        @test (
            sol.stats.nf, sol.stats.nsolve, sol.stats.nw, sol.stats.njacs,
            sol.stats.nnonliniter,
        ) == (41603, 1603, 800, 800, 1603)
        @test calls[] == sol.stats.nf
    end
end
