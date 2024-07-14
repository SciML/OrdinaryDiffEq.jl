using OrdinaryDiffEq, DiffEqDevTools, DiffEqBase, Test, Random
import ODEProblemLibrary: prob_ode_2Dlinear, prob_ode_linear

Random.seed!(100)
## Convergence Testing
dts = (1 / 2) .^ (7:-1:4)
testTol = 0.2

prob = prob_ode_linear
sol = solve(prob, OwrenZen3())
@test length(sol) < 20
sol = solve(prob, OwrenZen4())
@test length(sol) < 20
sol = solve(prob, OwrenZen5())
@test length(sol) < 20

sim = test_convergence(dts, prob, OwrenZen3(), dense_errors = true)
@test sim.𝒪est[:final]≈3 atol=testTol
@test sim.𝒪est[:L2]≈3 atol=testTol
sim = test_convergence(dts, prob, OwrenZen4(), dense_errors = true)
@test sim.𝒪est[:final]≈4 atol=testTol
@test sim.𝒪est[:L2]≈4 atol=testTol
sim = test_convergence(dts, prob, OwrenZen5(), dense_errors = true)
@test sim.𝒪est[:final]≈5 atol=testTol
@test sim.𝒪est[:L2]≈5 atol=testTol

prob = prob_ode_2Dlinear
sol = solve(prob, OwrenZen3())
@test length(sol) < 20
sol = solve(prob, OwrenZen4())
@test length(sol) < 20
sol = solve(prob, OwrenZen5())
@test length(sol) < 20

sim = test_convergence(dts, prob, OwrenZen3(), dense_errors = true)
@test sim.𝒪est[:final]≈3 atol=testTol
@test sim.𝒪est[:L2]≈3 atol=testTol
sim = test_convergence(dts, prob, OwrenZen4(), dense_errors = true)
@test sim.𝒪est[:final]≈4 atol=testTol
@test sim.𝒪est[:L2]≈4 atol=testTol
sim = test_convergence(dts, prob, OwrenZen5(), dense_errors = true)
@test sim.𝒪est[:final]≈5 atol=testTol
@test sim.𝒪est[:L2]≈5 atol=testTol
