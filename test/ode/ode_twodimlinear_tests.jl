using OrdinaryDiffEq
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_2Dlinear
prob = prob_ode_2Dlinear

integrator = init(prob,Euler();dt=1//2^(4))

integrator = init(prob,Tsit5();dt=1//2^(4))

solve!(integrator)

sol =solve(prob,Euler();dt=1//2^(4),maxiters=Inf)

sol =solve(prob,Tsit5();dt=1//2^(4))

# plot(sol,plot_analytic=true)

sol =solve(prob,ExplicitRK();dt=1//2^(4))
