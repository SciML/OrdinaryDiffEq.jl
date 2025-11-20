using OrdinaryDiffEq, ParameterizedFunctions, DiffEqDevTools, StaticArrays,
      OrdinaryDiffEqTaylorSeries
using Plots, BenchmarkTools, Random, LinearAlgebra
gr();

function lotka(du, u, p, t)
    x = u[1]
    y = u[2]
    du[1] = p[1] * x - p[2] * x * y
    du[2] = -p[3] * y + p[4] * x * y
end
prob_lotka = ODEProblem{true, SciMLBase.FullSpecialize}(
    lotka, [1.0, 1.0], (0.0, 10.0), [1.5, 1.0, 3.0, 1.0])
sol_lotka = solve(prob_lotka, Vern7(), abstol = 1e-14, reltol = 1e-14);

function fitzhugh(du, u, p, t)
    v = u[1]
    w = u[2]
    a = p[1]
    b = p[2]
    τinv = p[3]
    l = p[4]
    du[1] = v - v^3 / 3 - w + l
    du[2] = τinv * (v + a - b * w)
end
prob_fitzhugh = ODEProblem{true, SciMLBase.FullSpecialize}(
    fitzhugh, [1.0; 1.0], (0.0, 10.0), [0.7, 0.8, 1 / 12.5, 0.5])
sol_fitzhugh = solve(prob_fitzhugh, Vern7(), abstol = 1e-14, reltol = 1e-14);

function rigid(du, u, p, t)
    I₁ = p[1]
    I₂ = p[2]
    I₃ = p[3]
    du[1] = I₁ * u[2] * u[3]
    du[2] = I₂ * u[1] * u[3]
    du[3] = I₃ * u[1] * u[2] + 0.25 * sin(t)^2
end
prob_rigid = ODEProblem{true, SciMLBase.FullSpecialize}(
    rigid, [1.0; 0.0; 0.9], (0.0, 10.0), [-2.0, 1.25, -0.5])
sol_rigid = solve(prob_rigid, Vern7(), abstol = 1e-14, reltol = 1e-14);

function make_problem_random_dense(N; seed = 1)
    Random.seed!(seed)
    A = randn(N, N)
    f(du, u, p, t) = du .= A * u
    u0 = rand(N)
    ODEProblem{true, SciMLBase.FullSpecialize}(f, u0, (0.0, 10.0))
end

# 16 is fine, but Symbolics struggles with larger sizes like 64
prob_dense = make_problem_random_dense(16)
sol_dense = solve(prob_dense, Vern7(), abstol = 1e-14, reltol = 1e-14);

# Make work-precision plot
setups = [Dict(:alg => DP5())
          Dict(:alg => Tsit5())
          Dict(:alg => Vern6())
          Dict(:alg => Vern8())]
names = ["DP5", "Tsit5", "Vern6", "Vern8"]
for order in 6:2:12
    push!(names, "Taylor $(order)")
    push!(setups, Dict(:alg => ExplicitTaylor(order = Val(order + 1))))
end

function make_plot(
        prob, sol; numruns = 100, maxiters = 1000, absrange = (-13:-6), relrange = (-10:-3))
    abstols = 10.0 .^ absrange
    reltols = 10.0 .^ relrange
    wp = WorkPrecisionSet([prob], abstols, reltols, setups; names = names, appxsol = [sol],
        save_everystep = false, numruns = numruns, maxiters = maxiters)
    p = plot(wp)
    return p
end
p_lotka = make_plot(prob_lotka, sol_lotka)
p_fitzhugh = make_plot(prob_fitzhugh, sol_fitzhugh)
p_rigid = make_plot(prob_rigid, sol_rigid)
p_dense = make_plot(prob_dense, sol_dense)
