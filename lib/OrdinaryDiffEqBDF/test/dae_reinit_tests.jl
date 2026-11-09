using OrdinaryDiffEqBDF: DFBDF, DImplicitEuler
using SciMLBase: DAEProblem, init, reinit!, solve!, successful_retcode
using StaticArrays: SVector
using Test

@testset "DAE default and explicit-state reinitialization" begin
    function residual!(r, du, u, p, t)
        r[1] = du[1] + u[1]
        r[2] = u[2] - u[1]
        return nothing
    end
    residual(du, u, p, t) = [du[1] + u[1], u[2] - u[1]]
    static_residual(du, u, p, t) = SVector(du[1] + u[1], u[2] - u[1])
    scalar_residual(du, u, p, t) = du + u
    problems = (
        DAEProblem(residual!, [-1.0, -1.0], [1.0, 1.0], (0.0, 0.1)),
        DAEProblem(residual, [-1.0, -1.0], [1.0, 1.0], (0.0, 0.1)),
        DAEProblem(static_residual, SVector(-1.0, -1.0), SVector(1.0, 1.0), (0.0, 0.1)),
        DAEProblem(scalar_residual, -1.0, 1.0, (0.0, 0.1)),
    )
    for prob in problems, alg in (DImplicitEuler(), DFBDF()), reinit_dae in (true, false)
        @testset "$(typeof(prob.u0)), $(nameof(typeof(alg))), reinit_dae=$reinit_dae" begin
            u0, du0 = copy(prob.u0), copy(prob.du0)
            integrator = init(prob, alg; dt = 0.001, adaptive = false, abstol = 1.0e-8)
            sol = solve!(integrator)
            @test successful_retcode(sol)
            saved_initial = sol.u[1]
            for cycle in 1:2
                reinit!(integrator; reinit_dae)
                @test integrator.t == prob.tspan[1]
                @test integrator.u == u0
                @test integrator.sol.u[1] == u0
                sol = solve!(integrator)
                @test successful_retcode(sol)
                @test sol.t[end] == prob.tspan[2]
                @test maximum(abs.(sol.u[end] .- u0 .* exp(-sol.t[end]))) < 1.0e-4
                @test maximum(abs.(sol(0.0005) .- u0 .* exp(-0.0005))) < 1.0e-6
                @test saved_initial == u0
                @test prob.u0 == u0
                @test prob.du0 == du0
            end
            endpoint = copy(integrator.u)
            reinit!(integrator, endpoint; t0 = integrator.t, tf = 0.2)
            @test integrator.u == endpoint
            sol = solve!(integrator)
            @test successful_retcode(sol)
            @test sol.t[end] == 0.2
            @test sol.u[1] == endpoint
            @test maximum(abs.(sol.u[end] .- u0 .* exp(-sol.t[end]))) < 1.0e-4
            @test prob.u0 == u0
            @test prob.du0 == du0
        end
    end
end
