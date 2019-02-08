# This definitely needs cleaning
using OrdinaryDiffEq, DiffEqDevTools, Test, Random
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear
probArr = Vector{ODEProblem}(undef, 2)
probArr[1] = prob_ode_linear

probArr[2] = prob_ode_2Dlinear
Random.seed!(100)
## Convergence Testing
dts = 1 .//2 .^(8:-1:4)
testTol = 0.2



# Order Convergence test
for i = 1:2
  global dts
  prob = probArr[i]
  for j = 1:4
      sim = test_convergence(dts,prob,RichardsonEuler(j,j,j))
      @test sim.𝒪est[:final] ≈ j atol=testTol
   end
end
