using OrdinaryDiffEq,Plots,DiffEqProblemLibrary, DiffEqDevTools

prob = prob_ode_linear
println("Solve and Plot")
sol =solve(prob,Rosenbrock32())
TEST_PLOT && plot(sol,plot_analytic=true)
dt₀ = sol.t[2]

prob = prob_ode_2Dlinear

## Solve and plot
println("Solve and Plot")
tab = constructBogakiShampine3()
sol =solve(prob,ExplicitRK(),tableau=tab)
TEST_PLOT && plot(sol,plot_analytic=true)
dt₀ = sol.t[2]

bool1 = 1e-7 < dt₀ < .1
@test_throws ErrorException sol = solve(prob,Euler())
#TEST_PLOT && plot(sol,plot_analytic=true)
#dt₀ = sol.t[2]

tab = constructDormandPrince8_64bit()
sol3 =solve(prob,ExplicitRK(),tableau=tab)
TEST_PLOT && plot(sol3,plot_analytic=true)
dt₀ = sol3.t[2]

bool3 = 1e-7 < dt₀ < .3

bool1 && bool3
