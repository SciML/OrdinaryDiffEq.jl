using OrdinaryDiffEq, ParameterizedFunctions, DiffEqDevTools, StaticArrays,
      OrdinaryDiffEqTaylorSeries
using Plots, BenchmarkTools
gr();

abstols = 1.0 ./ 10.0 .^ (6:15)
reltols = 1.0 ./ 10.0 .^ (3:12);

function lotka(du, u, p, t)
    x = u[1]
    y = u[2]
    du[1] = p[1] * x - p[2] * x * y
    du[2] = -p[3] * y + p[4] * x * y
end
prob = ODEProblem{true, SciMLBase.FullSpecialize}(
    lotka, [1.0, 1.0], (0.0, 10.0), [1.5, 1.0, 3.0, 1.0])
sol = solve(prob, Vern7(), abstol = 1 / 10^14, reltol = 1 / 10^14)

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
wp = WorkPrecisionSet([prob], abstols, reltols, setups; names = names, appxsol = [sol],
    save_everystep = false, numruns = 100, maxiters = 10000)
plot(wp)
