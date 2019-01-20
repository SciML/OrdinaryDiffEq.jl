using OrdinaryDiffEq, DiffEqDevTools, Test, LinearAlgebra
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, van

testTol = 0.2

for prob in [prob_ode_linear]
  sim21 = test_convergence(1 .//2 .^(6:-1:3),prob,RadauIIA5())
  @test sim21.𝒪est[:final] ≈ 5 atol=testTol
end

# test adaptivity
vanstiff = ODEProblem{false}(van, [0;sqrt(3)], (0.0,1.0), 1e6)
@test length(solve(vanstiff, RadauIIA5())) < 110
@test length(solve(remake(vanstiff, p=1e7), RadauIIA5())) < 120
