using OrdinaryDiffEq

prob = prob_ode_linear
dt = 1//2^(4)
@time sol =solve(prob,cvode_BDF(),dt=dt,save_timeseries=true,)
@time sol =solve(prob,cvode_Adams(),dt=dt,save_timeseries=true)
@time sol =solve(prob,cvode_Adams(),dt=dt,save_timeseries=true,adaptive=false)

prob = prob_ode_2Dlinear
@time sol =solve(prob,cvode_BDF(),dt=dt,save_timeseries=true)
@time sol =solve(prob,cvode_Adams(),dt=dt,save_timeseries=true)
@time sol =solve(prob,cvode_Adams(),dt=dt,save_timeseries=true,adaptive=false)

length(sol)==17
