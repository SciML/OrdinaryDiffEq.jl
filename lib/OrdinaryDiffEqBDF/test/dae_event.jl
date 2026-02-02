using OrdinaryDiffEqBDF, Test

f = function (out, du, u, p, t)
    out[1] = -p[1] * u[1] + p[3] * u[2] * u[3] - du[1]
    out[2] = +p[1] * u[1] - p[2] * u[2]^2 - p[3] * u[2] * u[3] - du[2]
    return out[3] = u[1] + u[2] + u[3] - p[4]
end
u₀ = [1.0, 0, 0]
du₀ = [0.0, 0.0, 0.0]
p = [0.04, 3.0e7, 1.0e4, 1.0]
tspan = (0.0, 100.0)
differential_vars = [true, true, false]
prob = DAEProblem(f, du₀, u₀, tspan, p, differential_vars = differential_vars)
condition(u, t, integrator) = t in [50.0]
affect!(integrator) = integrator.p[4] = 2.0
cb = DiscreteCallback(condition, affect!)

#=
# Regression test
using Sundials
sol = solve(prob, IDA(), callback=cb, tstops=[50.0],abstol=1e-14,reltol=1e-14)
=#

p = [0.04, 3.0e7, 1.0e4, 1.0]
prob = DAEProblem(f, du₀, u₀, tspan, p, differential_vars = differential_vars)
sol = solve(prob, DFBDF(), callback = cb, tstops = [50.0], abstol = 1.0e-12, reltol = 1.0e-12)
@test sol.t[end] == 100.0
@test sol[end][1] ≈ 0.686300529575259 atol = 1.0e-7
@test sol[end][2] ≈ 2.0797982209353813e-6 atol = 1.0e-7
@test sol[end][3] ≈ 1.31369739062652 atol = 1.0e-7
