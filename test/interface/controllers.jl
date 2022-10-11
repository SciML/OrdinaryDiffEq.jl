using OrdinaryDiffEq, Test

f(u, p, t) = zero(u)
prob = ODEProblem(f, 1.0, (0.0, 10.0))

# Shouldn't error
sol = solve(prob, Tsit5(), controller = IController())
@test sol.retcode == :Success

sol = solve(prob, Tsit5(), controller = PIController(7 // 50, 2 // 25))
@test sol.retcode == :Success

sol = solve(prob, Tsit5(), controller = PIDController(0.7, -0.4))
@test sol.retcode == :Success

# OrdinaryDiffEq.jl#1703
# https://github.com/SciML/OrdinaryDiffEq.jl/issues/1703
prob = ODEProblem((du, u, p, t) -> du[1] = 1, [0.0], (0.0, 5.0))
sol = solve(prob, RDPK3SpFSAL49())
@test sol.retcode == :Success
