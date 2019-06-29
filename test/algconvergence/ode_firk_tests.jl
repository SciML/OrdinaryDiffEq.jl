using OrdinaryDiffEq, DiffEqDevTools, Test, LinearAlgebra
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, van

testTol = 0.2

for prob in [prob_ode_linear, prob_ode_2Dlinear]
  sim21 = test_convergence(1 .//2 .^(6:-1:3),prob,RadauIIA5())
  @test sim21.𝒪est[:final] ≈ 5 atol=testTol
end

# test adaptivity
for iip in (true, false)
  if iip
    vanstiff = ODEProblem{iip}(van, [0;sqrt(3)], (0.0,1.0), 1e6)
  else
    vanstiff = ODEProblem{false}((u,p,t)->van(u,p,t), [0;sqrt(3)], (0.0,1.0), 1e6)
  end
  sol = solve(vanstiff, RadauIIA5())
  if iip
    @test sol.destats.naccept + sol.destats.nreject > sol.destats.njacs # J reuse
    @test sol.destats.njacs < sol.destats.nw # W reuse
  end
  @test length(sol) < 150
  @test length(solve(remake(vanstiff, p=1e7), RadauIIA5())) < 150
  @test length(solve(remake(vanstiff, p=1e7), reltol=[1e-4, 1e-6], RadauIIA5())) < 170
  @test length(solve(remake(vanstiff, p=1e7), RadauIIA5(), reltol=1e-9, abstol=1e-9)) < 870
  @test length(solve(remake(vanstiff, p=1e9), RadauIIA5())) < 170
  @test length(solve(remake(vanstiff, p=1e10), RadauIIA5())) < 190
end
