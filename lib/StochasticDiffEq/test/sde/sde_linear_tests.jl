using StochasticDiffEq, DiffEqDevTools, Test, Random
using SDEProblemLibrary: prob_sde_linear
Random.seed!(100)
prob = prob_sde_linear

## Solve and plot
println("Solve and Plot")
sol = solve(prob, EM(), dt = 1 // 2^(4))
sol = solve(prob, RKMil(), dt = 1 // 2^(4))
sol = solve(prob, RKMilCommute(), dt = 1 // 2^(4))
sol = solve(prob, RKMilGeneral(), dt = 1 // 2^(4))
sol = solve(prob, SRI(), dt = 1 // 2^(4))
sol = solve(prob, SRIW1(), dt = 1 // 2^(4))
trajectories = 100
## Convergence Testing
println("Convergence Test on Linear")
dts = (1 // 2) .^ (9:-1:4) #14->7 good plot with higher num Monte

sim = test_convergence(dts, prob, EM(), trajectories = trajectories)

sim2 = test_convergence(dts, prob, RKMil(), trajectories = trajectories)

sim21 = test_convergence(dts, prob, RKMilGeneral(), trajectories = trajectories)

sim22 = test_convergence(dts, prob, RKMilCommute(), trajectories = trajectories)

sim3 = test_convergence(dts, prob, SRI(), trajectories = trajectories)

#TEST_PLOT && plot(plot(sim),plot(sim2),plot(sim3),layout=@layout([a b c]),size=(1200,600))

@test abs(sim.𝒪est[:l2] - 0.5) + abs(sim2.𝒪est[:l∞] - 1) + abs(sim3.𝒪est[:final] - 1.5) < 0.441  #High tolerance since low dts for testing!

@test abs(sim.𝒪est[:l2] - 0.5) + abs(sim21.𝒪est[:l∞] - 1) + abs(sim3.𝒪est[:final] - 1.5) < 0.441  #High tolerance since low dts for testing!

@test abs(sim22.𝒪est[:l∞] - 1) < 0.3

# test reinit
integrator = init(prob, EM(), dt = 1 // 2^(4))
solve!(integrator)
reinit!(integrator)
solve!(integrator)

# test reinit
prob2 = SDEProblem((u, p, t) -> prob.f(u, p, t), prob.g, prob.u0, prob.tspan)
integrator = init(prob2, EM(), dt = 1 // 2^(4), tstops = [1 // 2], saveat = [1 // 3])
solve!(integrator)
reinit!(integrator)
solve!(integrator)
