using StochasticDiffEq, Test

f1(u, p, t) = zero(u)
g1(u, p, t) = u isa Number ? one(u) : ones(size(u))
dt = 1 // 2^(4)
prob1 = SDEProblem{false}(f1, g1, 0.0, (0.0, 1.0))
integrator = init(prob1, EM(), dt = dt, save_noise = true)

k = integrator.W.dW
@test integrator.W.dW != 0
solve!(integrator)
@test integrator.sol.W[2] == k

prob1 = SDEProblem{false}(f1, g1, zeros(4), (0.0, 1.0))
sol = solve(prob1, EM(), dt = dt, save_noise = true)
@test sol.W[2] != zeros(4)

f1(du, u, p, t) = du .= 0.0
g1(du, u, p, t) = du .= 1.0
dt = 1 // 2^(4)
prob1 = SDEProblem(f1, g1, zeros(4), (0.0, 1.0))
sol = solve(prob1, EM(), dt = dt, save_noise = true)
@test sol.W[2] != zeros(4)
