using OrdinaryDiffEq,Plots
prob = prob_ode_2Dlinear
## Solve and plot
println("Solve and Plot")
sol =solve(prob,Rosenbrock32(),dt=1/2^4)

sol =solve(prob,ExplicitRK(tableau = constructBogakiShampine3()),dt=1/2^4)
val1 = maximum(abs.(sol.u[end] - sol.u_analytic[end]))
TEST_PLOT && plot(sol,plot_analytic=true)

sol2 =solve(prob,ExplicitRK(tableau = constructDormandPrince()),dt=1/2^4)
val2 = maximum(abs.(sol2.u[end] - sol2.u_analytic[end]))
TEST_PLOT && plot(sol2,plot_analytic=true)


sol3 =solve(prob,ExplicitRK(tableau = constructRKF8(Float64)),dt=1/2^4)
val3 = maximum(abs.(sol3.u[end] - sol3.u_analytic[end]))
TEST_PLOT && plot(sol3,plot_analytic=true)

@test length(sol.t)>length(sol2.t)>=length(sol3.t)
@test max(val1,val2,val3)<2e-3
