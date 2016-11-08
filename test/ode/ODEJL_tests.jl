using OrdinaryDiffEq

prob = prob_ode_linear
dt = 1/2^(4) #The initial timestepping size. It will automatically assigned if not given.

sol =solve(prob,feuler();dt=dt,save_timeseries=true)
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob,rk23(),dt=dt,save_timeseries=true)

sol =solve(prob,rk45(),dt=dt,save_timeseries=true)

sol =solve(prob,feh78(),dt=dt,save_timeseries=true)

sol =solve(prob,ModifiedRosenbrockIntegrator(),dt=dt,save_timeseries=true)

sol =solve(prob,midpoint(),dt=dt,save_timeseries=true)

sol =solve(prob,heun(),dt=dt,save_timeseries=true)

sol =solve(prob,rk4(),dt=dt,save_timeseries=true)

sol =solve(prob,feh45(),dt=dt,save_timeseries=true)

prob = prob_ode_2Dlinear

sol =solve(prob,feuler(),dt=dt,save_timeseries=true)
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob,rk23(),dt=dt,save_timeseries=true)

sol =solve(prob,rk45(),dt=dt,save_timeseries=true)

sol =solve(prob,feh78(),dt=dt,save_timeseries=true)

#sol =solve(prob,dt=0,save_timeseries=true,alg=:ode23s) #ODE.jl issues
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob,midpoint(),dt=dt,save_timeseries=true)

sol =solve(prob,heun(),dt=dt,save_timeseries=true)

sol =solve(prob,rk4(),dt=dt,save_timeseries=true)

sol =solve(prob,feh45(),dt=dt,save_timeseries=true)

#=
prob = prob_ode_bigfloat2Dlinear

sol =solve(prob,feuler(),dt=dt,save_timeseries=true)
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob,rk23(),dt=dt,save_timeseries=true)

sol =solve(prob,rk45(),dt=dt,save_timeseries=true)

sol =solve(prob,feh78(),dt=dt,save_timeseries=true)

#sol =solve(prob,dt=0,save_timeseries=true,alg=:ode23s) #ODE.jl issues
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob,midpoint(),dt=dt,save_timeseries=true)

sol =solve(prob,heun(),dt=dt,save_timeseries=true)

sol =solve(prob,rk4(),dt=dt,save_timeseries=true)

sol =solve(prob,feh45(),dt=dt,save_timeseries=true)
=#

true
