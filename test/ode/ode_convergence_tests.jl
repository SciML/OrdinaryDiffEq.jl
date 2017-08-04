# This definitely needs cleaning
using OrdinaryDiffEq, DiffEqDevTools, DiffEqBase,
      DiffEqProblemLibrary, Base.Test
probArr = Vector{ODEProblem}(2)
probArr[1] = prob_ode_linear

probArr[2] = prob_ode_2Dlinear
srand(100)
## Convergence Testing
println("Convergence Test on Linear")
dts = 1.//2.^(8:-1:4)
testTol = 0.2

for i = 1:2
  prob = probArr[i]
  println("Special RKs")
  sim = test_convergence(dts,prob,Euler())
  @test abs(sim.𝒪est[:final]-1) < testTol
  sim2 = test_convergence(dts,prob,Midpoint())
  @test abs(sim2.𝒪est[:l∞]-2) < testTol
  sim3 = test_convergence(dts,prob,RK4())
  @test abs(sim3.𝒪est[:l∞]-4) < testTol
  sim4 = test_convergence(dts,prob,BS3())
  @test abs(sim4.𝒪est[:l2]-3) < testTol

  ### Stiff Solvers

  println("Convergence Test on Stiff")
  dts = 1.//2.^(8:-1:4)

  sim11 = test_convergence(dts,prob,ImplicitEuler())
  @test abs(sim11.𝒪est[:final]-1) < testTol

  sim12 = test_convergence(dts,prob,
          GenericImplicitEuler(nlsolve=NLSOLVEJL_SETUP(autodiff=true)))
  @test abs(sim12.𝒪est[:final]-1) < testTol
  sim122 = test_convergence(dts,prob,
           GenericImplicitEuler(nlsolve=NLSOLVEJL_SETUP(autodiff=false)))

  sim13 = test_convergence(dts,prob,Trapezoid())
  @test abs(sim13.𝒪est[:final]-2) < testTol

  sim14 = test_convergence(dts,prob,
          GenericTrapezoid(nlsolve=NLSOLVEJL_SETUP(autodiff=true)))
  @test abs(sim14.𝒪est[:final]-2) < testTol
  sim142 = test_convergence(dts,prob,
           GenericTrapezoid(nlsolve=NLSOLVEJL_SETUP(autodiff=false)))
  @test abs(sim142.𝒪est[:final]-2) < testTol
end

prob = prob_ode_linear
dts = 1.//2.^(8:-1:4)
sim12 = test_convergence(dts,prob,ImplicitEuler(),dense_errors=true)
@test abs(sim12.𝒪est[:final]-1) < testTol
sim13 = test_convergence(dts,prob,Trapezoid(),dense_errors=true)
@test abs(sim13.𝒪est[:final]-2) < testTol

@time sol = solve(prob,ImplicitEuler(),dt=1/1000)
@time sol = solve(prob,GenericImplicitEuler(),dt=1/1000,adaptive=true)

@time sol = solve(prob,Trapezoid(),dt=1/1000)
@time sol = solve(prob,GenericTrapezoid(),dt=1/1000,adaptive=true)

prob = prob_ode_2Dlinear
dts = 1.//2.^(8:-1:4)
sim12 = test_convergence(dts,prob,ImplicitEuler())
@test abs(sim12.𝒪est[:final]-1) < testTol
sim13 = test_convergence(dts,prob,Trapezoid())
@test abs(sim13.𝒪est[:final]-2) < testTol

@time sol = solve(prob,ImplicitEuler(),dt=1/1000,reltol=1e-1)











@time sol = solve(prob,GenericImplicitEuler(),dt=1/1000,adaptive=true,reltol=1e-1)















@time sol = solve(prob,Trapezoid(),dt=1/1000)













@time sol = solve(prob,GenericTrapezoid(),dt=1/1000,adaptive=true,reltol=1e-3)
