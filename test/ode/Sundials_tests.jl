using OrdinaryDiffEq, DiffEqProblemLibrary

prob = prob_ode_linear
dt = 1//2^(4)
saveat = float(collect(0:dt:1))
@time sol = solve(prob,cvode_BDF)
@time sol = solve(prob,cvode_Adams)
@time sol = solve(prob,cvode_Adams,saveat=saveat)

prob = prob_ode_2Dlinear
@time sol = solve(prob,cvode_BDF)
@time sol = solve(prob,cvode_Adams)
@time sol = solve(prob,cvode_Adams,saveat=saveat)

length(sol)==17
true
