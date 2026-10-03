using OrdinaryDiffEqNordsieck
using OrdinaryDiffEqTsit5: Tsit5
using DiffEqDevTools, Test, LinearAlgebra
import ODEProblemLibrary: prob_ode_bigfloatlinear,
    prob_ode_bigfloat2Dlinear,
    prob_ode_linear, prob_ode_2Dlinear

probArr = [prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear]
testTol = 0.25
dts = 1 .// (2 .^ (10:-1:5))

@testset "Nordsieck Convergence Tests" begin
    for i in eachindex(probArr)
        sim = test_convergence(dts, probArr[i], AN5())
        @test sim.𝒪est[:final] ≈ 5 atol = testTol
        @test sim.𝒪est[:l2] ≈ 5 atol = testTol
        @test sim.𝒪est[:l∞] ≈ 5 atol = testTol
    end
end

@testset "JVODE corrector failure recovery" begin
    function lorenz!(du, u, p, t)
        du[1] = 10 * (u[2] - u[1])
        du[2] = u[1] * (28 - u[3]) - u[2]
        du[3] = u[1] * u[2] - (8 / 3) * u[3]
        return nothing
    end
    lorenz(u, p, t) = (du = similar(u); lorenz!(du, u, p, t); du)
    for f in (lorenz!, lorenz)
        prob = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 2.0))
        ref = solve(prob, Tsit5(); reltol = 1.0e-12, abstol = 1.0e-12)
        for alg in (JVODE_Adams(), JVODE_BDF())
            sol = solve(prob, alg; reltol = 1.0e-6, abstol = 1.0e-6)
            @test SciMLBase.successful_retcode(sol)
            @test sol.t[end] == 2.0
            @test sol.stats.nreject > 0
            @test norm(sol.u[end] - ref.u[end], Inf) < 1.0e-3
        end
    end
end

probArr = [
    prob_ode_linear,
    prob_ode_2Dlinear,
]
@testset "Nordsieck Adaptivity Tests: AN5" begin
    for i in eachindex(probArr)
        prob = probArr[i]
        sol = solve(prob, AN5(), reltol = 1.0e-6)
        @test length(sol.t) < 11
        @test SciMLBase.successful_retcode(sol)
        exact = prob.f.analytic(prob.u0, prob.p, prob.tspan[end])
        @test exact ≈ sol.u[end] atol = 1.0e-5
    end
end

@testset "Nordsieck Adaptivity Tests: JVODE" begin
    for i in eachindex(probArr),
            sol in [JVODE_Adams(), JVODE_BDF()]

        prob = probArr[i]
        sol = solve(prob, sol, reltol = 1.0e-4, abstol = 1.0e-7)
        @test length(sol.t) < 22
        @test SciMLBase.successful_retcode(sol)
        exact = prob.f.analytic(prob.u0, prob.p, prob.tspan[end])
        @test norm(exact - sol.u[end], Inf) < 3.0e-3
    end
end
