# This definitely needs cleaning
using OrdinaryDiffEq, DiffEqDevTools, DiffEqBase,
      DiffEqProblemLibrary, Base.Test
probArr = Vector{ODEProblem}(2)
probArr[1] = prob_ode_linear

probArr[2] = prob_ode_2Dlinear
srand(100)
## Convergence Testing
dts = 1.//2.^(8:-1:4)
dts1 = 1.//2.^(9:-1:5)
testTol = 0.2

for i = 1:2
  prob = probArr[i]
  sim = test_convergence(dts,prob,Euler())
  @test abs(sim.𝒪est[:final]-1) < testTol
  sim2 = test_convergence(dts,prob,Heun())
  @test abs(sim2.𝒪est[:l∞]-2) < testTol
  sim2 = test_convergence(dts,prob,Ralston())
  @test abs(sim2.𝒪est[:l∞]-2) < testTol
  sim2 = test_convergence(dts,prob,Midpoint())
  @test abs(sim2.𝒪est[:l∞]-2) < testTol
  sim3 = test_convergence(dts,prob,RK4())
  @test abs(sim3.𝒪est[:l∞]-4) < testTol
  sim4 = test_convergence(dts,prob,BS3())
  @test abs(sim4.𝒪est[:l2]-3) < testTol
  alg = CarpenterKennedy2N54()
  @test abs(test_convergence(dts,prob,alg).𝒪est[:l∞]-OrdinaryDiffEq.alg_order(alg)) < testTol
  sim5 = test_convergence(dts, prob, AB3())
  @test abs(sim5.𝒪est[:l2]-3) < testTol
  sim6 = test_convergence(dts,prob,ABM32())
  @test abs(sim6.𝒪est[:l2]-3) < testTol
  sim7 = test_convergence(dts, prob, AB4())
  @test abs(sim7.𝒪est[:l2]-4) < testTol
  sim8 = test_convergence(dts1,prob,ABM43())  #using dts1 due to floating point error in convergence test
  @test abs(sim8.𝒪est[:l2]-4) < testTol
  sim9 = test_convergence(dts,prob,AB5())
  @test abs(sim9.𝒪est[:l2]-5) < testTol
  sim10 = test_convergence(dts,prob,ABM54())
  @test abs(sim10.𝒪est[:l2]-5) < testTol
  sim101 = test_convergence(dts,prob,VCAB3())
  @test abs(sim101.𝒪est[:l2]-3) < testTol
  sim102 = test_convergence(dts,prob,VCAB4())
  @test abs(sim102.𝒪est[:l2]-4) < testTol
  sim103 = test_convergence(dts,prob,VCAB5())
  @test abs(sim103.𝒪est[:l2]-5) < testTol
  sim104 = test_convergence(dts,prob,VCABM3())
  @test abs(sim104.𝒪est[:l2]-3) < testTol
  sim105 = test_convergence(dts,prob,VCABM4())
  @test abs(sim105.𝒪est[:l2]-4) < testTol

  println("Stiff Solvers")

  dts = 1.//2.^(9:-1:5)

  sim11 = test_convergence(dts,prob,ImplicitEuler(extrapolant = :linear))
  @test abs(sim11.𝒪est[:final]-1) < testTol

  sim12 = test_convergence(dts,prob,
          GenericImplicitEuler(nlsolve=OrdinaryDiffEq.NLSOLVEJL_SETUP(autodiff=true)))
  @test abs(sim12.𝒪est[:final]-1) < testTol

  sim122 = test_convergence(dts,prob,
           GenericImplicitEuler(nlsolve=OrdinaryDiffEq.NLSOLVEJL_SETUP(autodiff=false)))

  sim13 = test_convergence(dts,prob,ImplicitMidpoint())
  @test abs(sim13.𝒪est[:final]-2) < testTol

  sim13 = test_convergence(dts,prob,Trapezoid())
  @test abs(sim13.𝒪est[:final]-2) < testTol

  sim14 = test_convergence(dts,prob,
          GenericTrapezoid(nlsolve=OrdinaryDiffEq.NLSOLVEJL_SETUP(autodiff=true)))
  @test abs(sim14.𝒪est[:final]-2) < testTol

  sim142 = test_convergence(dts,prob,
           GenericTrapezoid(nlsolve=OrdinaryDiffEq.NLSOLVEJL_SETUP(autodiff=false)))
  @test abs(sim142.𝒪est[:final]-2) < testTol

  sim14 = test_convergence(dts,prob,TRBDF2())
  @test abs(sim14.𝒪est[:final]-2) < testTol

  sim152 = test_convergence(dts,prob,TRBDF2(autodiff=false))
  @test abs(sim152.𝒪est[:final]-2) < testTol+0.1

  sim15 = test_convergence(dts,prob,SDIRK2())
  @test abs(sim15.𝒪est[:final]-2) < testTol

  sim152 = test_convergence(dts,prob,SSPSDIRK2())
  @test abs(sim152.𝒪est[:final]-2) < testTol

  sim16 = test_convergence(dts,prob,Kvaerno3())
  @test abs(sim16.𝒪est[:final]-3) < testTol

  sim17 = test_convergence(dts,prob,KenCarp3())
  @test abs(sim17.𝒪est[:final]-3) < testTol

  #####################################
  # BDF
  #####################################

  sim = test_convergence(dts,prob,ABDF2(extrapolant = :linear))
  @test abs(sim.𝒪est[:final]-2) < testTol
  @test abs(sim.𝒪est[:l2]-2) < testTol
  @test abs(sim.𝒪est[:l∞]-2) < testTol

  dts = 1.//2.^(7:-1:4)
  println("Higher Order")

  sim18 = test_convergence(dts,prob,Cash4())
  @test abs(sim18.𝒪est[:final]-4) < testTol

  sim19 = test_convergence(dts,prob,Hairer4())
  @test abs(sim19.𝒪est[:final]-4) < testTol

  sim110 = test_convergence(dts,prob,Hairer42())
  @test abs(sim110.𝒪est[:final]-4) < testTol

  sim111 = test_convergence(dts,prob,Kvaerno4())
  @test abs(sim111.𝒪est[:final]-4) < testTol

  sim112 = test_convergence(dts,prob,KenCarp4())
  @test abs(sim112.𝒪est[:final]-4) < testTol

  sim113 = test_convergence(dts,prob,Kvaerno5())
  @test abs(sim113.𝒪est[:final]-5) < testTol

  sim114 = test_convergence(dts,prob,KenCarp5())
  @test abs(sim114.𝒪est[:final]-5) < testTol
end
