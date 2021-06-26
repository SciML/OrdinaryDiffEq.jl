using OrdinaryDiffEq, DiffEqDevTools, Test, Random
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear

Random.seed!(100)

testTol = 0.25

f = (u,p,t)->cos(t)
prob_ode_sin = ODEProblem(ODEFunction(f; analytic=(u0,p,t)->sin(t)), 0.,(0.0,1.0))

f = (du,u,p,t)->du[1]=cos(t)
prob_ode_sin_inplace = ODEProblem(ODEFunction(f;analytic=(u0,p,t)->[sin(t)]), [0.], (0.0,1.0))

f = (u,p,t)->sin(u)
prob_ode_nonlinear = ODEProblem(ODEFunction(f;analytic=(u0,p,t)->2*acot(exp(-t)*cot(0.5))), 1.,(0.,0.5))

f = (du,u,p,t)->du[1]=sin(u[1])
prob_ode_nonlinear_inplace = ODEProblem(ODEFunction(f;analytic=(u0,p,t)->[2*acot(exp(-t)*cot(0.5))]),[1.],(0.,0.5))

test_problems_only_time = [prob_ode_sin, prob_ode_sin_inplace]
test_problems_linear = [prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear]
test_problems_nonlinear = [prob_ode_nonlinear, prob_ode_nonlinear_inplace]

# Test the memory usage, cf. #640
# Note: Basically, the size of the integrator should be the size of the cache
# plus the size of the initial condition (stored is integ.sol.prob.u0) if the
# keyword argument `alias_u0` is not set to `true` (default).
# Note: The memory requirements of the 2N methods can be reduced if an assignment
# of the form `tmp = A2end[i]*tmp + dt*f(u, p, t+c2end[i]*dt)` can be carried out
# without saving `f(u, p, t+c2end[i]*dt)` as `k`.
u0_large = rand(10^6)
prob_ode_large = ODEProblem((du,u,p,t)-> du .= u, u0_large, (0.0,1.0))

@testset "ORK256" begin
  alg = ORK256()
  alg2 = ORK256(;williamson_condition=false)
  dts = 1 ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 2
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CarpenterKennedy2N54" begin
  alg = CarpenterKennedy2N54()
  alg2 = CarpenterKennedy2N54(;williamson_condition=false)
  dts = 1 ./ 2 .^(7:-1:3)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 2
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "SHLDDRK64" begin
  alg = SHLDDRK64()
  alg2 = SHLDDRK64(;williamson_condition=true)
  dts = 1 ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 2
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "DGLDDRK73_C" begin
  alg = DGLDDRK73_C()
  alg2 = DGLDDRK73_C(;williamson_condition=false)
  dts = 1 ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 2
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "DGLDDRK84_C" begin
  alg = DGLDDRK84_C()
  alg2 = DGLDDRK84_C(;williamson_condition=false)
  dts = 1 ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 2
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "DGLDDRK84_F" begin
  alg = DGLDDRK84_F()
  alg2 = DGLDDRK84_F(;williamson_condition=false)
  dts = 1 ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 2
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "NDBLSRK124" begin
  alg = NDBLSRK124()
  alg2 = NDBLSRK124(;williamson_condition=false)
  dts = 1 ./ 2 .^(7:-1:3)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 2
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "NDBLSRK134" begin
  alg = NDBLSRK134()
  alg2 = NDBLSRK134(;williamson_condition=false)
  dts = 1 ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 2
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "NDBLSRK144" begin
  alg = NDBLSRK144()
  alg2 = NDBLSRK144(;williamson_condition=false)
  dts = 1 ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    sim = test_convergence(dts, prob, alg2)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 2
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  integ = init(prob_ode_large, alg2, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CFRLDDRK64" begin
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
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "TSLDDRK74" begin
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
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

# Methods from Carpenter, Kennedy, Lewis (2000)

function RemakeNew(p::ODEProblem)
  u1 = @. BigFloat(p.u0)
  tsp1 = @. BigFloat(p.tspan)
  remake(p; u0 = u1, tspan = tsp1)
end

test_problems_only_time_BigFloat = @. RemakeNew(test_problems_only_time)
test_problems_linear_BigFloat = @. RemakeNew(test_problems_linear)
f = (u,p,t)->sin(u)
prob_nonlinear_A = ODEProblem(ODEFunction(f;analytic=(u0,p,t)->2*acot(exp(-t)*cot(BigFloat(0.5)))), BigFloat(1.),(BigFloat(0.),BigFloat(0.5)))

f = (du,u,p,t)->du[1]=sin(u[1])
prob_nonlinear_B = ODEProblem(ODEFunction(f;analytic=(u0,p,t)->[2*acot(exp(-t)*cot(BigFloat(0.5)))]),[BigFloat(1.)],(BigFloat(0.),BigFloat(0.5)))
test_problems_nonlinear_BigFloat = [prob_nonlinear_A,prob_nonlinear_B]

@testset "CKLLSRK43_2" begin
  alg = CKLLSRK43_2()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+1 atol=testTol    # This scheme has linear order of 4
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 7
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK54_3C" begin
  alg = CKLLSRK54_3C()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ 1 atol=testTol          # The CI plot is linear but the evaluated order is 1
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 7
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK95_4S" begin
  alg = CKLLSRK95_4S()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 7
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK95_4C" begin
  alg = CKLLSRK95_4C()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 7
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK95_4M" begin
  alg = CKLLSRK95_4M()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 7
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK54_3C_3R" begin
  alg = CKLLSRK54_3C_3R()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 9
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 10
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 9
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK54_3M_3R" begin
  alg = CKLLSRK54_3M_3R()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+1 atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+0.5 atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 9
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 10
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 9
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK54_3N_3R" begin
  alg = CKLLSRK54_3N_3R()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 9
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 10
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 9
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK85_4C_3R" begin
  alg = CKLLSRK85_4C_3R()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 9
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 10
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 9
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK85_4M_3R" begin
  alg = CKLLSRK85_4M_3R()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 9
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 10
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 9
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK85_4P_3R" begin
  alg = CKLLSRK85_4P_3R()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+2 atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 9
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 10
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 9
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK54_3N_4R" begin
  alg = CKLLSRK54_3N_4R()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 11
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 12
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 11
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK54_3M_4R" begin
  alg = CKLLSRK54_3M_4R()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+0.5 atol=testTol                              # This scheme has linear orderof 4.5
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 11
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 12
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 11
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK65_4M_4R" begin
  alg = CKLLSRK65_4M_4R()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 11
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 12
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 11
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK85_4FM_4R" begin
  alg = CKLLSRK85_4FM_4R()
  dts = BigFloat(1) ./ 2 .^(10:-1:6)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) + 1 atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 11
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 12
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 11
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "CKLLSRK75_4M_5R" begin
  alg = CKLLSRK75_4M_5R()
  dts = BigFloat(1) ./ 2 .^(8:-1:4)
  for prob in test_problems_only_time_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_linear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear_BigFloat
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, adaptive=false,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 13
  integ = init(prob_ode_large, alg, adaptive=true,dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 14
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 13
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end


# Methods from Parsani, Ketcheson, Deconinck (2013)

@testset "ParsaniKetchesonDeconinck3S32" begin
  alg = ParsaniKetchesonDeconinck3S32()
  dts = 1 ./ 2 .^(7:-1:3)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    # higher order as pure quadrature
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+1 atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "ParsaniKetchesonDeconinck3S82" begin
  alg = ParsaniKetchesonDeconinck3S82()
  dts = 1 ./ 2 .^(8:-1:5)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    # higher order as pure quadrature
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+1 atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "ParsaniKetchesonDeconinck3S53" begin
  alg = ParsaniKetchesonDeconinck3S53()
  dts = 1 ./ 2 .^(7:-1:3)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    # higher order as pure quadrature
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+1 atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "ParsaniKetchesonDeconinck3S173" begin
  alg = ParsaniKetchesonDeconinck3S173()
  dts = 1 ./ 2 .^(7:-1:3)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    # higher order as pure quadrature
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+1 atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  dts = 1 ./ 2 .^(6:-1:3)
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=1
  end

  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "ParsaniKetchesonDeconinck3S94" begin
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
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "ParsaniKetchesonDeconinck3S184" begin
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
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "ParsaniKetchesonDeconinck3S105" begin
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
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "ParsaniKetchesonDeconinck3S205" begin
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
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

# Methods from Ranocha, Dalcin, Parsani, Ketcheson (2021)

@testset "RDPK3Sp35" begin
  alg = RDPK3Sp35()
  dts = 1 ./ 2 .^(7:-1:3)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    # higher order as pure quadrature
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+1 atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "RDPK3Sp49" begin
  alg = RDPK3Sp49()
  dts = 1 ./ 2 .^(5:-1:2)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    # higher order as pure quadrature
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+1 atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  dts = 1 ./ 2 .^(8:-1:2)
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "RDPK3Sp510" begin
  alg = RDPK3Sp510()
  dts = 1 ./ 2 .^(4.5:-1:1.5)
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
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "RDPK3SpFSAL35" begin
  alg = RDPK3SpFSAL35()
  dts = 1 ./ 2 .^(7:-1:3)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    # higher order as pure quadrature
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+1 atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "RDPK3SpFSAL49" begin
  alg = RDPK3SpFSAL49()
  dts = 1 ./ 2 .^(5:-1:2)
  for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    # higher order as pure quadrature
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+1 atol=testTol
  end
  for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  dts = 1 ./ 2 .^(8:-1:2)
  for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
  end
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end

@testset "RDPK3SpFSAL510" begin
  alg = RDPK3SpFSAL510()
  dts = 1 ./ 2 .^(4.5:-1:1.5)
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
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
  integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
  @test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
  # test whether aliasing u0 is bad
  new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
  sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
  sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
  @test sol_old[end] ≈ sol_new[end]
end
