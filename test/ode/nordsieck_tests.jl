using OrdinaryDiffEq, DiffEqDevTools, Base.Test

srand(100)
linear_bigαN = big"0.5"
f_linearbig = (u,p,t) -> (linear_bigαN*u)
f_2dlinearbig = (du,u,p,t) -> (du.=linear_bigαN*u)
(f::typeof(f_linearbig))(::Type{Val{:analytic}},u0,p,t) = u0*exp(linear_bigαN*t)
(f::typeof(f_2dlinearbig))(::Type{Val{:analytic}},u0,p,t) = u0*exp.(linear_bigαN*t)
probArr = [ODEProblem(f_2dlinearbig, big.(rand(4,2)), (0,1.)),
           ODEProblem(f_linearbig, big"0.5", (0,1.))]
testTol = 0.2
dts = 1.//2.^(10:-1:6)

for prob in probArr
  println("AN5")
  sim = test_convergence(dts,prob,AN5())
  @test abs(sim.𝒪est[:final]-5) < testTol
  @test abs(sim.𝒪est[:l2]-5) < testTol
  @test abs(sim.𝒪est[:l∞]-5) < testTol
end
