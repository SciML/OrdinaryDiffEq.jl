# Introduction to the ODE Solvers

using OrdinaryDiffEq

prob = prob_ode_linear
dt = 1//2^(4) #The initial timestepping size. It will automatically assigned if not given.
tspan = [0,1] # The timespan. This is the default if not given.
sol =solve(prob,dopri5(),dt=dt,save_timeseries=true)
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob,dop853();dt=1//2^(4),save_timeseries=true)

sol =solve(prob,odex();dt=1//2^(4),save_timeseries=true)

sol =solve(prob,seulex();dt=1//2^(4),save_timeseries=true)

sol =solve(prob,radau();dt=1//2^(4),save_timeseries=true)

sol =solve(prob,radau5();dt=1//2^(4),save_timeseries=true)

prob = prob_ode_2Dlinear

sol =solve(prob,dopri5(),dt=dt,save_timeseries=true)

sol =solve(prob,dop853();dt=1//2^(4),save_timeseries=true)

sol =solve(prob,odex();dt=1//2^(4),save_timeseries=true)

sol =solve(prob,seulex();dt=1//2^(4),save_timeseries=true)

sol =solve(prob,radau();dt=1//2^(4),save_timeseries=true)

sol =solve(prob,radau5();dt=1//2^(4),save_timeseries=true)

#=
prob = prob_ode_bigfloat2Dlinear

sol =solve(prob,dopri5(),dt=dt,save_timeseries=true)

sol =solve(prob,dop853();dt=1//2^(4),save_timeseries=true)

sol =solve(prob,odex();dt=1//2^(4),save_timeseries=true)

sol =solve(prob,seulex();dt=1//2^(4),save_timeseries=true)

sol =solve(prob,radau();dt=1//2^(4),save_timeseries=true)

sol =solve(prob,radau5();dt=1//2^(4),save_timeseries=true)
=#

true
