# Introduction to the ODE Solvers

using OrdinaryDiffEq

prob = prob_ode_linear
dt = 1//2^(4) #The initial timestepping size. It will automatically assigned if not given.
tspan = [0;1.0] # The timespan. This is the default if not given.
sol =solve(prob,dopri5,dt=dt)
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob,dop853;dt=1//2^(4))

sol =solve(prob,odex;dt=1//2^(4))

sol =solve(prob,seulex;dt=1//2^(4))

sol =solve(prob,radau;dt=1//2^(4))

sol =solve(prob,radau5;dt=1//2^(4))

prob = prob_ode_2Dlinear

sol =solve(prob,dopri5,dt=dt)

sol =solve(prob,dop853;dt=1//2^(4))

sol =solve(prob,odex;dt=1//2^(4))

sol =solve(prob,seulex;dt=1//2^(4))

sol =solve(prob,radau;dt=1//2^(4))

sol =solve(prob,radau5;dt=1//2^(4))

#=
prob = prob_ode_bigfloat2Dlinear

sol =solve(prob,dopri5,dt=dt)

sol =solve(prob,dop853;dt=1//2^(4))

sol =solve(prob,odex;dt=1//2^(4))

sol =solve(prob,seulex;dt=1//2^(4))

sol =solve(prob,radau;dt=1//2^(4))

sol =solve(prob,radau5;dt=1//2^(4))
=#

true
