using OrdinaryDiffEq, DiffEqDevTools, Test, Random, Plots
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear
using Quadmath

Random.seed!(100)

dts = Float128.(1 .// 2 .^ (6:-1:2))
testTol = 0.35

# Tests on simple problem

f = (u, p, t) -> cos(t)
prob_ode_sin = ODEProblem(
                    ODEFunction(f; analytic = (u0, p, t) -> sin(t)), 
                    Float128(0.0), 
                    (Float128(0.0), Float128(1.0)))

f = (du, u, p, t) -> du[1] = cos(t)
prob_ode_sin_inplace = ODEProblem(
                    ODEFunction(f; analytic = (u0, p, t) -> [sin(t)]), 
                    [Float128(0.0)], 
                    (Float128(0.0), Float128(1.0)))

f = (u, p, t) -> sin(u)
prob_ode_nonlinear = ODEProblem(
                    ODEFunction(f; analytic = (u0, p, t) -> Float128(2.0) * acot(exp(-t) 
                    * cot(Float128(0.5)))),
                    Float128(1.0),
                    (Float128(0.0), Float128(0.5)))

f = (du, u, p, t) -> du[1] = sin(u[1])
prob_ode_nonlinear_inplace = ODEProblem(
                    ODEFunction(f; analytic = (u0, p, t) -> [Float128(2.0) * acot(exp(-t) 
                    * cot(Float128(0.5)))]),
                    [Float128(1.0)], 
                    (Float128(0.0), Float128(0.5)))

test_problems_only_time = [prob_ode_sin, prob_ode_sin_inplace]
test_problems_linear = [prob_ode_bigfloat2Dlinear]
test_problems_nonlinear = [prob_ode_nonlinear, prob_ode_nonlinear_inplace]


println("QPRK98")
alg = QPRK98()

for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    sim.𝒪est[:final]
    @test sim.𝒪est[:final]≈OrdinaryDiffEq.alg_order(alg)+1 atol=testTol
    sol = solve(prob, alg, adaptive = true, save_everystep = false)
    sol_exact = prob.f.analytic(prob.u0, prob.p, sol.t[end])
    @test minimum(abs.(sol.u[end] .- sol_exact) .< 1e-12)
end

for prob in test_problems_linear
    sim = test_convergence(BigFloat.(dts), prob, alg)
    sim.𝒪est[:final]
    @test sim.𝒪est[:final]≈OrdinaryDiffEq.alg_order(alg)+1 atol=testTol
    sol = solve(prob, alg, adaptive = true, save_everystep = false)
    sol_exact = prob.f.analytic(prob.u0, prob.p, sol.t[end])
    @test minimum(abs.(sol.u[end] .- sol_exact) .< 1e-8)
end

for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    sim.𝒪est[:final]
    @test sim.𝒪est[:final]≈OrdinaryDiffEq.alg_order(alg)+2.5 atol=testTol
    sol = solve(prob, alg, adaptive = true, save_everystep = false)
    sol_exact = prob.f.analytic(prob.u0, prob.p, sol.t[end])
    @test minimum(abs.(sol.u[end] .- sol_exact) .< 1e-11)
end


