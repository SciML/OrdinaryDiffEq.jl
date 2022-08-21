using OrdinaryDiffEq, DiffEqBase, Test, DiffEqDevTools,
      Random

import ODEProblemLibrary: prob_ode_bigfloatlinear,
                                               prob_ode_bigfloat2Dlinear,
                                               prob_ode_2Dlinear

## Convergence Testing
println("Convergence Test on Linear")

testTol = 1
prob = prob_ode_2Dlinear
println("Feagin RKs")
dts = (1//2) .^ (4:-1:2)
sol = solve(prob,Feagin10(),dt=dts[1])
prob = remake(prob_ode_bigfloat2Dlinear,tspan=(big(0)//1,big(1)//1))
sol = solve(prob,Feagin10(),dt=dts[1])

prob = remake(prob_ode_bigfloat2Dlinear,tspan=(big(0.0),big(1.0)))
dts = (1//2) .^ (4:-1:2)
sim = test_convergence(dts,prob,Feagin10())
@test abs(sim.𝒪est[:final]-8) < testTol #Lowered due to low test dt

sim = test_convergence(dts,prob,Feagin12())
@test abs(sim.𝒪est[:final]-12) < testTol

sim = test_convergence(dts,prob,Feagin14())
@test abs(sim.𝒪est[:final]-15) < testTol #Upped to 15 for test

prob = prob_ode_bigfloatlinear

dts = (1//2) .^ (6:-1:3)
sim = test_convergence(dts,prob,Feagin10())
@test abs(sim.𝒪est[:final]-10) < testTol

dts = (1//2) .^ (4:-1:2)
sim = test_convergence(dts,prob,Feagin12())
@test abs(sim.𝒪est[:final]-12) < testTol

sim = test_convergence(dts,prob,Feagin14())
@test abs(sim.𝒪est[:final]-15) < testTol #Upped to 15 for test

prob = prob_ode_bigfloat2Dlinear

#compile
sol =solve(prob,Feagin10(),dt=dts[1])
sol =solve(prob,Feagin12(),dt=dts[1])
sol =solve(prob,Feagin14(),dt=dts[1])

#test
@time sol =solve(prob,Feagin10(),dt=dts[1])
@time sol =solve(prob,Feagin12(),dt=dts[1])
@time sol =solve(prob,Feagin14(),dt=dts[1])
