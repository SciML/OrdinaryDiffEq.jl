# Introduction to the ODE Solvers

using OrdinaryDiffEq

prob = prob_ode_linear
dt = 1//2^(4) #The initial timestepping size. It will automatically assigned if not given.
tspan = [0,1] # The timespan. This is the default if not given.
sol =solve(prob::ODEProblem,dopri5(),dt=dt,save_timeseries=true)
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob::ODEProblem,dop853();dt=1//2^(4),save_timeseries=true)

sol =solve(prob::ODEProblem,odex();dt=1//2^(4),save_timeseries=true)

sol =solve(prob::ODEProblem,seulex();dt=1//2^(4),save_timeseries=true)

sol =solve(prob::ODEProblem,radau();dt=1//2^(4),save_timeseries=true)

sol =solve(prob::ODEProblem,radau5();dt=1//2^(4),save_timeseries=true)

prob = prob_ode_2Dlinear

sol =solve(prob::ODEProblem,dopri5(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,dop853();dt=1//2^(4),save_timeseries=true)

sol =solve(prob::ODEProblem,odex();dt=1//2^(4),save_timeseries=true)

sol =solve(prob::ODEProblem,seulex();dt=1//2^(4),save_timeseries=true)

sol =solve(prob::ODEProblem,radau();dt=1//2^(4),save_timeseries=true)

sol =solve(prob::ODEProblem,radau5();dt=1//2^(4),save_timeseries=true)

#=
prob = prob_ode_bigfloat2Dlinear

sol =solve(prob::ODEProblem,dopri5(),dt=dt,save_timeseries=true)

sol =solve(prob::ODEProblem,dop853();dt=1//2^(4),save_timeseries=true)

sol =solve(prob::ODEProblem,odex();dt=1//2^(4),save_timeseries=true)

sol =solve(prob::ODEProblem,seulex();dt=1//2^(4),save_timeseries=true)

sol =solve(prob::ODEProblem,radau();dt=1//2^(4),save_timeseries=true)

sol =solve(prob::ODEProblem,radau5();dt=1//2^(4),save_timeseries=true)
=#

true
