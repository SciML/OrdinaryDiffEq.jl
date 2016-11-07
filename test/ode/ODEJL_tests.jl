using OrdinaryDiffEq

prob = prob_ode_linear
dt = 1/2^(4) #The initial timestepping size. It will automatically assigned if not given.

sol =solve(prob::ODEProblem,feuler();dt=dt,save_timeseries=true)
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob::ODEProblem,rk23(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,rk45(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,feh78(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,ModifiedRosenbrockIntegrator(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,midpoint(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,heun(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,rk4(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,feh45(),dt=dt,save_timeseries=true)

prob = prob_ode_2Dlinear

sol =solve(prob::ODEProblem,feuler(),dt=dt,save_timeseries=true)
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob::ODEProblem,rk23(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,rk45(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,feh78(),dt=dt,save_timeseries=true)

#sol =solve(prob::ODEProblem,dt=0,save_timeseries=true,alg=:ode23s) #ODE.jl issues
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob::ODEProblem,midpoint(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,heun(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,rk4(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,feh45(),dt=dt,save_timeseries=true)

#=
prob = prob_ode_bigfloat2Dlinear

sol =solve(prob::ODEProblem,feuler(),dt=dt,save_timeseries=true)
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob::ODEProblem,rk23(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,rk45(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,feh78(),dt=dt,save_timeseries=true)

#sol =solve(prob::ODEProblem,dt=0,save_timeseries=true,alg=:ode23s) #ODE.jl issues
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob::ODEProblem,midpoint(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,heun(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,rk4(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,feh45(),dt=dt,save_timeseries=true)
=#

true
