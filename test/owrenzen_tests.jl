using OrdinaryDiffEq, DiffEqDevTools, DiffEqBase,
      DiffEqProblemLibrary, Base.Test
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

sim = test_convergence(dts,prob,OwrenZen3())
@test abs(sim.𝒪est[:final]-3) < testTol
sim = test_convergence(dts,prob,OwrenZen4())
@test abs(sim.𝒪est[:final]-4) < testTol
sim = test_convergence(dts,prob,OwrenZen5())
@test abs(sim.𝒪est[:final]-5) < testTol

prob = prob_ode_2Dlinear
sol = solve(prob,OwrenZen3())
@test length(sol) < 20
sol = solve(prob,OwrenZen4())
@test length(sol) < 20
sol = solve(prob,OwrenZen5())
@test length(sol) < 20

sim = test_convergence(dts,prob,OwrenZen3())
@test abs(sim.𝒪est[:final]-3) < testTol
sim = test_convergence(dts,prob,OwrenZen4())
@test abs(sim.𝒪est[:final]-4) < testTol
sim = test_convergence(dts,prob,OwrenZen5())
@test abs(sim.𝒪est[:final]-5) < testTol
