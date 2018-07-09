## Breakout these since no other test of their adaptivity

using OrdinaryDiffEq, DiffEqDevTools, Test
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear,
       prob_ode_bigfloat2Dlinear

dts = 1.//2.^(6:-1:3)
testTol = 0.2

const linear_bigα3 = parse(BigFloat,"1.01")
f_2dlinearbig = (du,u,p,t) -> begin
  for i in 1:length(u)
    du[i] = linear_bigα3*u[i]
  end
end
(f::typeof(f_2dlinearbig))(::Type{Val{:analytic}},u0,p,t) = u0*exp.(1.01*t)
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

prob = prob_ode_bigfloat2Dlinear

sim = test_convergence(dts,prob,Rosenbrock23(linsolve=LinSolveFactorize(qrfact!)))
@test abs(sim.𝒪est[:final]-2) < testTol

sol = solve(prob,Rosenbrock23(linsolve=LinSolveFactorize(qrfact!)))
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

### ROS3P()

prob = prob_ode_linear

sim = test_convergence(dts,prob,ROS3P())
@test abs(sim.𝒪est[:final]-3) < testTol

sol = solve(prob,ROS3P())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,ROS3P())
@test abs(sim.𝒪est[:final]-3) < testTol

sol = solve(prob,ROS3P())
@test length(sol) < 20

### Rodas3()

prob = prob_ode_linear

sim = test_convergence(dts,prob,Rodas3())
@test abs(sim.𝒪est[:final]-3) < testTol

sol = solve(prob,Rodas3())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Rodas3())
@test abs(sim.𝒪est[:final]-3) < testTol

sol = solve(prob,Rodas3())
@test length(sol) < 20

println("4th order Rosenbrocks")

### RosShamp4

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

### Veldd4

prob = prob_ode_linear

sim = test_convergence(dts,prob,Veldd4())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,Veldd4())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Veldd4())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,Veldd4())
@test length(sol) < 20

### Velds4

prob = prob_ode_linear

sim = test_convergence(dts,prob,Velds4())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,Velds4())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Velds4())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,Velds4())
@test length(sol) < 20

### GRK4T

prob = prob_ode_linear

sim = test_convergence(dts,prob,GRK4T())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,GRK4T())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,GRK4T())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,GRK4T())
@test length(sol) < 20

### GRK4A
dts = 1.//2.^(7:-1:4)

prob = prob_ode_linear

sim = test_convergence(dts,prob,GRK4A())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,GRK4A())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,GRK4A())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,GRK4A())
@test length(sol) < 20

### Ros4LStab

prob = prob_ode_linear

sim = test_convergence(dts,prob,Ros4LStab())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,Ros4LStab())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Ros4LStab())
@test abs(sim.𝒪est[:final]-4) < testTol

sol = solve(prob,Ros4LStab())
@test length(sol) < 20

### Rodas4 Algorithms

println("RODAS")

dts = 1.//2.^(7:-1:4)

prob = prob_ode_linear

sim = test_convergence(dts,prob,Rodas4(),dense_errors=true)
@test abs(sim.𝒪est[:final]-4) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol

sol = solve(prob,Rodas4())
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas4(autodiff=false),dense_errors=true)
@test abs(sim.𝒪est[:final]-4) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol

sol = solve(prob,Rodas4(autodiff=false))
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas42(),dense_errors=true)
@test abs(sim.𝒪est[:final]-5) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol

sol = solve(prob,Rodas42())
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas4P(),dense_errors=true)
@test abs(sim.𝒪est[:final]-4) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol

sol = solve(prob,Rodas4P())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Rodas4(),dense_errors=true)
@test abs(sim.𝒪est[:final]-4) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol

sol = solve(prob,Rodas4())
@test length(sol) < 20

println("Rodas4 with finite diff")

sim = test_convergence(dts,prob,Rodas4(autodiff=false),dense_errors=true)
@test abs(sim.𝒪est[:final]-4) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol

sol = solve(prob,Rodas4(autodiff=false))
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas4(autodiff=false,
                       diff_type=Val{:forward}),dense_errors=true)
@test abs(sim.𝒪est[:final]-4) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol

sol = solve(prob,Rodas4(autodiff=false,diff_type=Val{:forward}))
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas4(autodiff=false,
                       diff_type=Val{:complex}),dense_errors=true)
@test abs(sim.𝒪est[:final]-4) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol

sol = solve(prob,Rodas4(autodiff=false,diff_type=Val{:complex}))
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas42(),dense_errors=true)
@test abs(sim.𝒪est[:final]-5) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol

sol = solve(prob,Rodas42())
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas4P(),dense_errors=true)
@test abs(sim.𝒪est[:final]-4) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol

sol = solve(prob,Rodas4P())
@test length(sol) < 20

### Rodas5
println("Rodas5")

prob = prob_ode_linear

dts = 1.//2.^(7:-1:3)
sim = test_convergence(dts,prob,Rodas5(),dense_errors=true)
@test abs(sim.𝒪est[:final]-5) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol

sol = solve(prob,Rodas5())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Rodas5(),dense_errors=true)
@test abs(sim.𝒪est[:final]-5) < testTol
@test abs(sim.𝒪est[:L2]-4) < testTol

sol = solve(prob,Rodas5())
@test length(sol) < 20
