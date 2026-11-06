using OrdinaryDiffEqFeagin, DiffEqBase, Test, DiffEqDevTools,
    Random, LinearAlgebra
using JLArrays

import ODEProblemLibrary: prob_ode_bigfloatlinear,
    prob_ode_bigfloat2Dlinear,
    prob_ode_2Dlinear

@testset "alg_adaptive_order extends OrdinaryDiffEqCore" begin
    @test OrdinaryDiffEqFeagin.OrdinaryDiffEqCore.alg_adaptive_order(Feagin10()) == 8
    @test OrdinaryDiffEqFeagin.OrdinaryDiffEqCore.alg_adaptive_order(Feagin14()) == 12
end

## Convergence Testing
println("Convergence Test on Linear")

testTol = 1
prob = prob_ode_2Dlinear
println("Feagin RKs")
dts = (1 // 2) .^ (4:-1:2)
sol = solve(prob, Feagin10(), dt = dts[1])
prob = remake(prob_ode_bigfloat2Dlinear, tspan = (big(0) // 1, big(1) // 1))
sol = solve(prob, Feagin10(), dt = dts[1])

prob = remake(prob_ode_bigfloat2Dlinear, tspan = (big(0.0), big(1.0)))
dts = (1 // 2) .^ (4:-1:2)
sim = test_convergence(dts, prob, Feagin10())
@test abs(sim.𝒪est[:final] - 8) < testTol #Lowered due to low test dt

sim = test_convergence(dts, prob, Feagin12())
@test abs(sim.𝒪est[:final] - 12) < testTol

sim = test_convergence(dts, prob, Feagin14())
@test abs(sim.𝒪est[:final] - 15) < testTol #Upped to 15 for test

prob = prob_ode_bigfloatlinear

dts = (1 // 2) .^ (6:-1:3)
sim = test_convergence(dts, prob, Feagin10())
@test abs(sim.𝒪est[:final] - 10) < testTol

dts = (1 // 2) .^ (4:-1:2)
sim = test_convergence(dts, prob, Feagin12())
@test abs(sim.𝒪est[:final] - 12) < testTol

sim = test_convergence(dts, prob, Feagin14())
@test abs(sim.𝒪est[:final] - 15) < testTol #Upped to 15 for test

prob = prob_ode_bigfloat2Dlinear

#compile
sol = solve(prob, Feagin10(), dt = dts[1])
sol = solve(prob, Feagin12(), dt = dts[1])
sol = solve(prob, Feagin14(), dt = dts[1])

#test
@time sol = solve(prob, Feagin10(), dt = dts[1])
@time sol = solve(prob, Feagin12(), dt = dts[1])
@time sol = solve(prob, Feagin14(), dt = dts[1])

@testset "Feagin in-place JLArray under allowscalar(false)" begin
    Ahost = [-1.0 0.2 0.0; -0.1 -0.8 0.1; 0.0 -0.2 -0.5]
    A = JLArray(Ahost)
    function linear_iip!(du, u, A, t)
        mul!(du, A, u)
        return nothing
    end
    initial = [1.0, 0.0, 0.0]
    u0 = JLArray(initial)
    prob_jl = ODEProblem(linear_iip!, u0, (0.0, 1.0), A)
    reference = exp(Ahost) * initial
    JLArrays.allowscalar(false)
    for alg in (Feagin10(), Feagin12(), Feagin14()), adaptive in (true, false)
        sol = solve(prob_jl, alg; adaptive, dt = 0.01, abstol = 1.0e-8, reltol = 1.0e-8)
        @test successful_retcode(sol)
        @test eltype(sol.u[end]) === Float64
        @test sol.u[end] isa JLArray
        @test Array(sol.u[end]) ≈ reference atol = 1.0e-7 rtol = 1.0e-7
    end
end
