using OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test

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
