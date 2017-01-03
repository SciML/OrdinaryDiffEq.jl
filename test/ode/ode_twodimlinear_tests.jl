using OrdinaryDiffEq, DiffEqProblemLibrary
prob = prob_ode_2Dlinear

sol =solve(prob,Euler();dt=1//2^(4))

integrator = init(prob,Euler();dt=1//2^(4))

# plot(sol,plot_analytic=true)

sol =solve(prob,ExplicitRK();dt=1//2^(4))
