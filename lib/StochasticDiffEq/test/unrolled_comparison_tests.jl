using StochasticDiffEq, DiffEqDevTools, Test, Random
using SDEProblemLibrary: prob_sde_linear, prob_sde_2Dlinear
prob = prob_sde_linear

Random.seed!(100)
sol1 = solve(prob, SRI(), dt = 1 // 2^(4))
Random.seed!(100)
sol2 = solve(prob, SRIW1(), dt = 1 // 2^(4))

@test sol1.u[end] ≈ sol2.u[end]

prob = prob_sde_2Dlinear

Random.seed!(100)
sol1 = solve(prob, SRI(), dt = 1 // 2^(4))
Random.seed!(100)
sol2 = solve(prob, SRIW1(), dt = 1 // 2^(4))

@test sol1.u[end] ≈ sol2.u[end]

prob = prob_sde_linear

Random.seed!(100)
integrator1 = init(prob, SRI(), dt = 1 // 2^(4))
step!(integrator1);
step!(integrator1)

Random.seed!(100)
integrator2 = init(prob, SRIW1(), dt = 1 // 2^(4))
step!(integrator2);
step!(integrator2)

@test integrator1.EEst ≈ integrator2.EEst

prob = prob_sde_2Dlinear

Random.seed!(100)
integrator1 = init(prob, SRI(), dt = 1 // 2^(4))
step!(integrator1);
step!(integrator1)

Random.seed!(100)
integrator2 = init(prob, SRIW1(), dt = 1 // 2^(4))
step!(integrator2);
step!(integrator2)

@test integrator1.EEst ≈ integrator2.EEst
