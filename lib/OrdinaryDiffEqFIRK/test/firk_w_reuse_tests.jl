using OrdinaryDiffEqFIRK, Test
using SciMLBase: successful_retcode

@testset "RadauIIA3 iteration matrix reuse" begin
    function rotation!(du, u, p, t)
        ω = p * sum(abs2, u)
        du[1] = ω * u[2]
        du[2] = -ω * u[1]
        return nothing
    end

    frequency = 4π * 1.0e12
    prob = ODEProblem(rotation!, [1.0, 0.0], (0.0, 10.3e-12), frequency)
    sol = solve(prob, RadauIIA3(); reltol = 1.0e-6, abstol = 1.0e-9, dt = 1.0e-14)
    @test successful_retcode(sol)
    @test sol.stats.nw - sol.stats.njacs < sol.stats.naccept ÷ 10
end
