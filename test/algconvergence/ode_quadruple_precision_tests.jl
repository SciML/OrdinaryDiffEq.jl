using OrdinaryDiffEq, DiffEqDevTools, Test, Random
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear

Random.seed!(100)

dts = 1 .// 2 .^ (8:-1:4)
testTol = 0.25

f = (u, p, t) -> cos(t)
prob_ode_sin = ODEProblem(ODEFunction(f; analytic = (u0, p, t) -> sin(t)), 0.0, (0.0, 1.0))

f = (du, u, p, t) -> du[1] = cos(t)
prob_ode_sin_inplace = ODEProblem(ODEFunction(f; analytic = (u0, p, t) -> [sin(t)]), [0.0],
    (0.0, 1.0))

f = (u, p, t) -> sin(u)
prob_ode_nonlinear = ODEProblem(
    ODEFunction(f;
        analytic = (u0, p, t) -> 2 * acot(exp(-t) *
                                          cot(0.5))), 1.0,
    (0.0, 0.5))

f = (du, u, p, t) -> du[1] = sin(u[1])
prob_ode_nonlinear_inplace = ODEProblem(
    ODEFunction(f;
        analytic = (u0, p, t) -> [
            2 * acot(exp(-t) * cot(0.5))
        ]),
    [1.0], (0.0, 0.5))

test_problems_only_time = [prob_ode_sin, prob_ode_sin_inplace]
test_problems_linear = [prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear]
test_problems_nonlinear = [prob_ode_nonlinear, prob_ode_nonlinear_inplace]

f_ssp = (u, p, t) -> begin
    sin(10t) * u * (1 - u)
end
test_problem_ssp = ODEProblem(f_ssp, 0.1, (0.0, 8.0))
test_problem_ssp_long = ODEProblem(f_ssp, 0.1, (0.0, 1.e3))

f_ssp_inplace = (du, u, p, t) -> begin
    @. du = sin(10t) * u * (1 - u)
end
test_problem_ssp_inplace = ODEProblem(f_ssp_inplace, rand(3, 3), (0.0, 8.0))


println("QPRK98")
alg = QPRK98()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final]≈OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final]≈OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final]≈OrdinaryDiffEq.alg_order(alg) atol=testTol
end

#=
# test SSP coefficient
sol = solve(test_problem_ssp_long, alg, dt = OrdinaryDiffEq.ssp_coefficient(alg),
    dense = false)
@test all(sol.u .>= 0)
# test SSP property of dense output
sol = solve(test_problem_ssp, alg, dt = 1.0)
@test mapreduce(t -> all(0 .<= sol(t) .<= 1), (u, v) -> u && v,
    range(0, stop = 8, length = 50), init = true)
sol = solve(test_problem_ssp_inplace, alg, dt = 1.0)
@test mapreduce(t -> all(0 .<= sol(t) .<= 1), (u, v) -> u && v,
    range(0, stop = 8, length = 50), init = true)
# test storage
integ = init(prob_ode_large, alg, dt = 1.e-2, save_start = false, save_end = false,
    save_everystep = false)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
integ = init(prob_ode_large, alg, dt = 1.e-2, save_start = false, save_end = false,
    save_everystep = false, alias_u0 = true)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
=#