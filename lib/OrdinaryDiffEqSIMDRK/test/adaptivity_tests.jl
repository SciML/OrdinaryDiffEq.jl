using OrdinaryDiffEqSIMDRK, StaticArrays, Test

function lorenz(u, p, t)
    return SA[10.0(u[2] - u[1]), u[1] * (28.0 - u[3]) - u[2], u[1] * u[2] - (8 / 3) * u[3]]
end
u0 = SA[1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz, u0, tspan)
sol = solve(prob, MER5v2())
@test length(sol.t) < 1100

sol = solve(prob, MER6v2())
@test length(sol.t) < 1100

sol = solve(prob, RK6v4())
@test length(sol.t) < 1700
