using OrdinaryDiffEq

prob = prob_ode_linear
sol =solve(prob,dopri5(),dt=1//2^(4))
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob,dop853();dt=1//2^(4))

sol =solve(prob,odex();dt=1//2^(4))

sol =solve(prob,seulex();dt=1//2^(4))

sol =solve(prob,radau();dt=1//2^(4))

sol =solve(prob,radau5();dt=1//2^(4))

prob = prob_ode_2Dlinear

sol =solve(prob,dopri5(),dt=dt)

sol =solve(prob,dop853();dt=1//2^(4))

sol =solve(prob,odex();dt=1//2^(4))

sol =solve(prob,seulex();dt=1//2^(4))

sol =solve(prob,radau();dt=1//2^(4))

sol =solve(prob,radau5();dt=1//2^(4))

#=
prob = prob_ode_bigfloat2Dlinear

sol =solve(prob,dopri5(),dt=dt)

sol =solve(prob,dop853();dt=1//2^(4))

sol =solve(prob,odex();dt=1//2^(4))

sol =solve(prob,seulex();dt=1//2^(4))

sol =solve(prob,radau();dt=1//2^(4))

sol =solve(prob,radau5();dt=1//2^(4))
=#
