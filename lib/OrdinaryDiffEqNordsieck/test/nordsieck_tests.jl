using OrdinaryDiffEq, DiffEqDevTools, Test, LinearAlgebra
import ODEProblemLibrary: prob_ode_bigfloatlinear,
                          prob_ode_bigfloat2Dlinear,
                          prob_ode_linear, prob_ode_2Dlinear

probArr = [prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear]
testTol = 0.25
dts = 1 .// (2 .^ (10:-1:5))

@testset "Nordsieck Convergence Tests" begin
    for i in eachindex(probArr)
        sim = test_convergence(dts, probArr[i], AN5())
        @test sim.𝒪est[:final]≈5 atol=testTol
        @test sim.𝒪est[:l2]≈5 atol=testTol
        @test sim.𝒪est[:l∞]≈5 atol=testTol
    end
end

probArr = [prob_ode_linear,
    prob_ode_2Dlinear]
@testset "Nordsieck Adaptivity Tests: AN5" begin
    for i in eachindex(probArr)
        prob = probArr[i]
        sol = solve(prob, AN5(), reltol = 1e-6)
        @test length(sol.t) < 11
        exact = prob.f.analytic(prob.u0, prob.p, prob.tspan[end])
        @test exact≈sol[end] atol=1e-5
    end
end

@testset "Nordsieck Adaptivity Tests: JVODE" begin
    for i in eachindex(probArr),
        sol in [JVODE_Adams(), JVODE_BDF()]

        prob = probArr[i]
        sol = solve(prob, sol, reltol = 1e-4, abstol = 1e-7)
        @test length(sol.t) < 22
        exact = prob.f.analytic(prob.u0, prob.p, prob.tspan[end])
        @test norm(exact - sol[end], Inf) < 3e-3
    end
end