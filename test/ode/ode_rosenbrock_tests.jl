## Breakout these since no other test of their adaptivity

using OrdinaryDiffEq, DiffEqProblemLibrary, DiffEqDevTools, Base.Test

dts = 1.//2.^(8:-1:4)
testTol = 0.2

const linear_bigα = parse(BigFloat,"1.01")
f_2dlinearbig = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = linear_bigα*u[i]
  end
end
(f::typeof(f_2dlinearbig))(::Type{Val{:analytic}},t,u0) = u0*exp.(1.01*t)
prob_ode_bigfloat2Dlinear = ODEProblem(f_2dlinearbig,map(BigFloat,rand(4,2)).*ones(4,2)/2,(0.0,1.0))

### Rosenbrock23()

prob = prob_ode_linear

sim = test_convergence(dts,prob,Rosenbrock23())
@test abs(sim.𝒪est[:final]-2) < testTol

sol = solve(prob,Rosenbrock23())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Rosenbrock23())
@test abs(sim.𝒪est[:final]-2) < testTol

sol = solve(prob,Rosenbrock23())
@test length(sol) < 20

sim = test_convergence(dts,prob,Rosenbrock23())
@test abs(sim.𝒪est[:final]-2) < testTol

sol = solve(prob,Rosenbrock23())
@test length(sol) < 20



### Rosenbrock32()

prob = prob_ode_linear

sim = test_convergence(dts,prob,Rosenbrock32())
@test abs(sim.𝒪est[:final]-3) < testTol

sol = solve(prob,Rosenbrock32())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Rosenbrock32())
@test abs(sim.𝒪est[:final]-3) < testTol

sol = solve(prob,Rosenbrock32())
@test length(sol) < 20

prob = prob_ode_bigfloat2Dlinear

sim = test_convergence(dts,prob,Rosenbrock32(linsolve=LinSolveFactorize(qrfact!)))
@test abs(sim.𝒪est[:final]-3) < testTol

sol = solve(prob,Rosenbrock32())
@test length(sol) < 20

### RosShamp4

dts = 1.//2.^(8:-1:3)

prob = prob_ode_linear

sim = test_convergence(dts,prob,RosShamp4())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,RosShamp4())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,RosShamp4())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,RosShamp4())
@test length(sol) < 20

prob = prob_ode_bigfloat2Dlinear

sim = test_convergence(dts,prob,RosShamp4())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,RosShamp4())
@test length(sol) < 20

### Veldd4

dts = 1.//2.^(8:-1:3)

prob = prob_ode_linear

sim = test_convergence(dts,prob,Veldd4())
@test abs(sim.𝒪est[:final]-2) < testTol

sol = solve(prob,Veldd4())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Veldd4())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,Veldd4())
@test length(sol) < 20

prob = prob_ode_bigfloat2Dlinear

sim = test_convergence(dts,prob,Veldd4())
@test abs(sim.𝒪est[:final]-5.87) < testTol

sol = solve(prob,Veldd4())
@test length(sol) < 20

### Velds4

dts = 1.//2.^(8:-1:3)

prob = prob_ode_linear

sim = test_convergence(dts,prob,Velds4())
@test abs(sim.𝒪est[:final]-2) < testTol

sol = solve(prob,Velds4())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Velds4())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,Velds4())
@test length(sol) < 20

prob = prob_ode_bigfloat2Dlinear

sim = test_convergence(dts,prob,Velds4())
@test abs(sim.𝒪est[:final]-5.87) < testTol

sol = solve(prob,Velds4())
@test length(sol) < 20

### Test on Stiff

prob = deepcopy(prob_ode_rober)
prob.tspan = (0.0,1e5)

sol = solve(prob,Rosenbrock23())
