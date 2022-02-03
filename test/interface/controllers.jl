using OrdinaryDiffEq, Test

f(u,p,t) = zero(u)
prob = ODEProblem(f, 1.0, (0.0, 10.0))

# Shouldn't error
sol = solve(prob, Tsit5(), controller=IController())
@test sol.retcode === :Success

sol = solve(prob, Tsit5(), controller=PIController(7//50, 2//25))
@test sol.retcode === :Success

sol = solve(prob, Tsit5(), controller=PIDController(0.7, -0.4))
@test sol.retcode === :Success
