using OrdinaryDiffEq, Test
f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0.0, Inf)
prob = ODEProblem(f, u0, tspan)

# Shouldn't error, but should go unstable and abort
sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8, adaptive = false, dt = 0.1)
@test sol.retcode == :Unstable
