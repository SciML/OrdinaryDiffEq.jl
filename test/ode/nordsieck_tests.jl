using OrdinaryDiffEq, DiffEqDevTools, DiffEqProblemLibrary, Base.Test

probArr = [prob_ode_bigfloatlinear,
           prob_ode_bigfloat2Dlinear]
testTol = 0.2
dts = 1.//2.^(10:-1:6)

for prob in probArr
  println("AN5")
  sim = test_convergence(dts,prob,AN5())
  @test abs(sim.𝒪est[:final]-5) < testTol
  @test abs(sim.𝒪est[:l2]-5) < testTol
  @test abs(sim.𝒪est[:l∞]-5) < testTol
end
