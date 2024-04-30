using OrdinaryDiffEq, DiffEqDevTools,  Test

import ODEProblemLibrary: prob_ode_linear,
                          prob_ode_2Dlinear


probnum = prob_ode_linear
prob = prob_ode_2Dlinear

dts = (1 / 2) .^ (7:-1:4)
testTol = 0.2

### Tsit5()

println("Tsit5")
dts = (1 / 2) .^ (7:-1:3)
sim = test_convergence(dts, probnum, Tsit5())
@test abs.(sim.𝒪est[:l2] - 5) < testTol + 0.2
sim = test_convergence(dts, prob, Tsit5())
@test abs.(sim.𝒪est[:l2] - 5) < testTol + 0.2

### Tsit5 with relaxation

println("Tsit5 with relaxation")
dts = (1 / 2) .^ (7:-1:3)
sim = test_convergence(dts, probnum, Tsit5_for_relaxation())
@test abs.(sim.𝒪est[:l2] - 5) < testTol + 0.2
#sim = test_convergence(dts, prob, Tsit5_for_relaxation())
#@test abs.(sim.𝒪est[:l2] - 5) < testTol + 0.2