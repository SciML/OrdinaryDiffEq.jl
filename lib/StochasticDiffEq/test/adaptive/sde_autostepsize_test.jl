using StochasticDiffEq, Random
using SDEProblemLibrary: prob_sde_2Dlinear
Random.seed!(100)

#Let the solver determine the initial stepsize for you!
sol = solve(prob_sde_2Dlinear, SRI())

# using Plots
# plot(sol,plot_analytic=true)
# gui()

#Make sure it does a good job
sol.t[2] > 1.0e-7

sol = solve(prob_sde_2Dlinear, DRI1())
#Make sure it does a good job
sol.t[2] > 1.0e-7
