## Breakout these since no other test of their adaptivity

using OrdinaryDiffEq, DiffEqProblemLibrary, DiffEqDevTools

dts = 1.//2.^(8:-1:4)
testTol = 0.2
bools = Vector{Bool}(0)

### Rosenbrock23

prob = prob_ode_linear

sim = test_convergence(dts,prob,Rosenbrock23)
push!(bools,abs(sim.𝒪est[:final]-2) < testTol)

sol = solve(prob,Rosenbrock23)
push!(bools,length(sol) < 20)

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Rosenbrock23)
push!(bools,abs(sim.𝒪est[:final]-2) < testTol)

sol = solve(prob,Rosenbrock23)
push!(bools,length(sol) < 20)

### Rosenbrock32

prob = prob_ode_linear

sim = test_convergence(dts,prob,Rosenbrock32)
push!(bools,abs(sim.𝒪est[:final]-3) < testTol)

sol = solve(prob,Rosenbrock32)
push!(bools,length(sol) < 20)

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Rosenbrock32)
push!(bools,abs(sim.𝒪est[:final]-3) < testTol)

sol = solve(prob,Rosenbrock32)
push!(bools,length(sol) < 20)

### Test on Stiff

prob = deepcopy(prob_ode_rober)
prob.tspan = (0.0,1e5)

sol = solve(prob,Rosenbrock23)


minimum(bools)
