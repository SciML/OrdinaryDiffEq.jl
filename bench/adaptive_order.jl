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

# Profile
cache1 = init(prob, ExplicitTaylorAdaptiveOrder(min_order = Val(6), max_order = Val(12)))
OrdinaryDiffEqTaylorSeries.perform_step!(cache1, cache1.cache, false)
@profview for _ in 1:10000000
    OrdinaryDiffEqTaylorSeries.perform_step!(cache1, cache1.cache, false)
end
@btime OrdinaryDiffEqTaylorSeries.perform_step!(cache1, cache1.cache, false)

# Make work-precision plot
setups = [Dict(:alg => DP5())
          Dict(:alg => Tsit5())
          Dict(:alg => Vern6())
          Dict(:alg => Vern8())
          Dict(:alg => ExplicitTaylorAdaptiveOrder(min_order = Val(8), max_order = Val(12)))];
names = ["DP5", "Tsit5", "Vern6", "Vern8", "Taylor"]
for order in 6:2:12
    push!(names, "Taylor $(order)")
    push!(setups, Dict(:alg => ExplicitTaylor(order = Val(order + 1))))
end
wp = WorkPrecisionSet([prob], abstols, reltols, setups; names = names, appxsol = [sol],
    save_everystep = false, numruns = 100, maxiters = 10000)
plot(wp)
