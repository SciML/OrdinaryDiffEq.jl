using OrdinaryDiffEq, DiffEqProblemLibrary, DiffEqCallbacks, Base.Test

prob = prob_ode_2Dlinear

integrator = init(prob,Tsit5();dt=1//2^(4))
dt = integrator.dt
solve!(integrator)

u = copy(integrator.sol.u)
t = copy(integrator.sol.t)

reinit!(integrator)
integrator.dt = dt
solve!(integrator)

@test u == integrator.sol.u
@test t == integrator.sol.t

# Implicitly test if dt resets
integrator = init(prob,Tsit5())
dt = integrator.dt
solve!(integrator)

u = copy(integrator.sol.u)
t = copy(integrator.sol.t)

reinit!(integrator)
integrator.dt
solve!(integrator)

@test u == integrator.sol.u
@test t == integrator.sol.t

integrator = init(prob,Tsit5();dt=1//2^(4),tstops = [1//2], saveat = [1//4])
dt = integrator.dt
solve!(integrator)

u = copy(integrator.sol.u)
t = copy(integrator.sol.t)

reinit!(integrator)
integrator.dt = dt
solve!(integrator)

@test u == integrator.sol.u
@test t == integrator.sol.t

#callback test
g(t,u) = 2.0*t-2.0
u0=0.0
tspan = (0.0,4.0)
prob = ODEProblem(g,u0,tspan)
saved_values = SavedValues(Float64, Float64)
cb = SavingCallback((t,u,integrator)->u, saved_values)
integrator = init(prob,Tsit5();dt=1//2^(4),callback = cb)
dt = integrator.dt
solve!(integrator)

u = saved_values.saveval
t = saved_values.t
resize!(saved_values.t,0)
resize!(saved_values.saveval,0)
reinit!(integrator)
integrator.dt = dt
solve!(integrator)

@test u == saved_values.saveval
@test t == saved_values.t