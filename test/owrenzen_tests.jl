using OrdinaryDiffEq, DiffEqDevTools, DiffEqBase, Test, Random
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_2Dlinear, prob_ode_linear

srand(100)
## Convergence Testing
dts = 1.//2.^(7:-1:4)
testTol = 0.2



prob = prob_ode_linear
sol = solve(prob,OwrenZen3())
@test length(sol) < 20
sol = solve(prob,OwrenZen4())
@test length(sol) < 20
sol = solve(prob,OwrenZen5())
@test length(sol) < 20

sim = test_convergence(dts,prob,OwrenZen3(),dense_errors=true)
@test abs(sim.𝒪est[:final]-3) < testTol
@test abs(sim.𝒪est[:L2]-3) < testTol
sim = test_convergence(dts,prob,OwrenZen4(),dense_errors=true)
@test abs(sim.𝒪est[:final]-4) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol
sim = test_convergence(dts,prob,OwrenZen5(),dense_errors=true)
@test abs(sim.𝒪est[:final]-5) < testTol
@test abs(sim.𝒪est[:L2]-5) < testTol

prob = prob_ode_2Dlinear
sol = solve(prob,OwrenZen3())
@test length(sol) < 20
sol = solve(prob,OwrenZen4())
@test length(sol) < 20
sol = solve(prob,OwrenZen5())
@test length(sol) < 20

sim = test_convergence(dts,prob,OwrenZen3(),dense_errors=true)
@test abs(sim.𝒪est[:final]-3) < testTol
@test abs(sim.𝒪est[:L2]-3) < testTol
sim = test_convergence(dts,prob,OwrenZen4(),dense_errors=true)
@test abs(sim.𝒪est[:final]-4) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol
sim = test_convergence(dts,prob,OwrenZen5(),dense_errors=true)
@test abs(sim.𝒪est[:final]-5) < testTol
@test abs(sim.𝒪est[:L2]-5) < testTol
