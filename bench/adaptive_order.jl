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

function f(du, u, p, t)
    A = 20
    σ = 0.02
    t0 = 5.0
    q = 1 + A * exp(-((t - t0)^2) / (2 * σ^2))
    du[1] = - q * u[1]
    du[2] = - q * u[2]
end

prob2 = ODEProblem{true, SciMLBase.FullSpecialize}(f, [100.0, 100.0], (0.0, 10.0), nothing)
sol2 = solve(prob2, Vern7(), abstol = 1 / 10^14, reltol = 1 / 10^14);

@btime sol = solve(prob2, ExplicitTaylorAdaptiveOrder(min_order = Val(6), max_order = Val(12)), abstol = 1e-12, reltol = 1e-12);
cache = init(prob2, ExplicitTaylorAdaptiveOrder(min_order = Val(6), max_order = Val(12)), abstol = 1e-12, reltol = 1e-12);
sol = solve!(cache)

@btime sol = solve(prob2, ExplicitTaylor(order = Val(12)), abstol = 1e-12, reltol = 1e-12);
@btime sol = solve(prob2, ExplicitTaylor(order = Val(11)), abstol = 1e-12, reltol = 1e-12);

solve!(cache)
# Profile
function profile1()
    cache1 = init(prob, ExplicitTaylor(order = Val(8)))
    OrdinaryDiffEqTaylorSeries.perform_step!(cache1, cache1.cache, false)
    @profview for _ in 1:10000000
        OrdinaryDiffEqTaylorSeries.perform_step!(cache1, cache1.cache, false)
    end
    @btime OrdinaryDiffEqTaylorSeries.perform_step!($cache1, $cache1.cache, false)
end

profile1()

# Make work-precision plot"DP5", "Tsit5", "Vern6", "Vern8",
# setups = [Dict(:alg => DP5())
#           Dict(:alg => Tsit5())
#           Dict(:alg => Vern6())
#           Dict(:alg => Vern8())
setups = Any[Dict(:alg => ExplicitTaylorAdaptiveOrder(min_order = Val(6), max_order = Val(12)))];
names = ["Taylor"]
for order in 6:2:12
    push!(names, "Taylor $(order)")
    push!(setups, Dict(:alg => ExplicitTaylor(order = Val(order + 1))))
end
wp = WorkPrecisionSet([prob2], abstols, reltols, setups; names = names, appxsol = [sol2],
    save_everystep = false, numruns = 100, maxiters = 10000)
plot(wp)
