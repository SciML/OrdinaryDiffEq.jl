using OrdinaryDiffEq, DiffEqDevTools, Test, Random
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear

Random.seed!(100)

testTol = 0.25

f = (u,p,t)->cos(t)
(::typeof(f))(::Type{Val{:analytic}},u0,p,t) = sin(t)
prob_ode_sin = ODEProblem(f, 0.,(0.0,1.0))

f = (du,u,p,t)->du[1]=cos(t)
(::typeof(f))(::Type{Val{:analytic}},u0,p,t) = [sin(t)]
prob_ode_sin_inplace = ODEProblem(f, [0.], (0.0,1.0))

f = (u,p,t)->sin(u)
(::typeof(f))(::Type{Val{:analytic}},u0,p,t) = 2*acot(exp(-t)*cot(0.5))
prob_ode_nonlinear = ODEProblem(f, 1.,(0.,0.5))

f = (du,u,p,t)->du[1]=sin(u[1])
(::typeof(f))(::Type{Val{:analytic}},u0,p,t) = [2*acot(exp(-t)*cot(0.5))]
prob_ode_nonlinear_inplace = ODEProblem(f,[1.],(0.,0.5))

test_problems_only_time = [prob_ode_sin, prob_ode_sin_inplace]
test_problems_linear = [prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear]
test_problems_nonlinear = [prob_ode_nonlinear, prob_ode_nonlinear_inplace]


alg = ORK256()
dts = 1 ./ 2 .^(8:-1:4)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end


alg = CarpenterKennedy2N54()
dts = 1 ./ 2 .^(7:-1:3)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end


alg = LDDRK64()
dts = 1 ./ 2 .^(8:-1:4)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end


alg = CFRLDDRK64()
dts = 1 ./ 2 .^(7:-1:4)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end


alg = TSLDDRK74()
dts = 1 ./ 2 .^(8:-1:4)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end


alg = DGLDDRK73_C()
dts = 1 ./ 2 .^(8:-1:4)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end


alg = DGLDDRK84_C()
dts = 1 ./ 2 .^(8:-1:4)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end



alg = NDBLSRK124()
dts = 1 ./ 2 .^(7:-1:3)
for prob in test_problems_only_time
	sim = test_convergence(dts, prob, alg)
	@test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
	sim = test_convergence(dts, prob, alg)
	@test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
	sim = test_convergence(dts, prob, alg)
	@test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end


alg = NDBLSRK134()
dts = 1 ./ 2 .^(8:-1:4)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end


alg = NDBLSRK144()
dts = 1 ./ 2 .^(8:-1:4)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end


alg = ParsaniKetchesonDeconinck3S94()
dts = 1 ./ 2 .^(7:-1:3)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end


alg = ParsaniKetchesonDeconinck3S184()
dts = 1 ./ 2 .^(6:-1:2)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
dts = 1 ./ 2 .^(7:-1:2)
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end


alg = ParsaniKetchesonDeconinck3S105()
dts = 1 ./ 1.95 .^(5:-1:1)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
dts = 1 ./ 2 .^(5:-1:2)
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
dts = 1.5 ./ 2 .^(5:-1:2)
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end


alg = ParsaniKetchesonDeconinck3S205()
dts = 1 ./ 1.95 .^(5:-1:1)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
dts = 1 ./ 2 .^(5:-1:2)
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
dts = 1.5 ./ 2 .^(5:-1:2)
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
