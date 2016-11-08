using OrdinaryDiffEq

prob = prob_ode_linear
dt = 1//2^(4)
@time sol =solve(prob,cvode_BDF,dt=dt)
@time sol =solve(prob,cvode_Adams,dt=dt)
@time sol =solve(prob,cvode_Adams,dt=dt,adaptive=false)

prob = prob_ode_2Dlinear
@time sol =solve(prob,cvode_BDF,dt=dt)
@time sol =solve(prob,cvode_Adams,dt=dt)
@time sol =solve(prob,cvode_Adams,dt=dt,adaptive=false)

length(sol)==17
