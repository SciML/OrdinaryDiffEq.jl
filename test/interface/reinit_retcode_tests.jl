using OrdinaryDiffEq, Test
f(u, p, t) = 2 * u
u_positive(u, t, integrator) = u > 0
terminate_if_u_pos = DiscreteCallback(u_positive, terminate!)

prob = ODEProblem(f, 1.0, (0.0, 1.0)) # positive initial condition > positive u > :Terminated
integrator = init(prob, Tsit5(), callback = terminate_if_u_pos)
sol1 = solve!(integrator)
@test sol1.retcode == :Terminated

reinit!(integrator, -1.0) # negative initial condition > negative u > :Success!
sol2 = solve!(integrator)
@test sol2.retcode == :Success

prob = ODEProblem(f, -1.0, (0.0, 1.0)) # negative initial condition > negative u > :Success!
integrator = init(prob, Tsit5(), callback = terminate_if_u_pos)
sol3 = solve!(integrator)
@test sol3.retcode == :Success

# https://github.com/SciML/DifferentialEquations.jl/issues/904
f(u, p, t) = -u
prob = ODEProblem(f, 1.0, (0.0, 1.0), save_on = false, save_start = false)
int = init(prob, Tsit5())

reinit!(int, 0.0)
solve!(int)
