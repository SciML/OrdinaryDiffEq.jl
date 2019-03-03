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

# Test the memory usage, cf. #640
# Note: Basically, the size of the integrator should be the size of the cache
# plus the size of the initial condition (stored is integ.sol.prob.u0) if the
# keyword argument `alias_u0` is not set to `true` (default).
# Note: The memory requirements of the 2N methods can be reduced if an assignment
# of the form `tmp = A2end[i]*tmp + dt*f(u, p, t+c2end[i]*dt)` can be carried out
# without saving `f(u, p, t+c2end[i]*dt)` as `k`.
u0_large = rand(10^6)
prob_ode_large = ODEProblem((du,u,p,t)-> du .= u, u0_large, (0.0,1.0))


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
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
# test whether aliasing u0 is bad
new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
@test sol_old[end] ≈ sol_new[end]


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
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
# test whether aliasing u0 is bad
new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
@test sol_old[end] ≈ sol_new[end]


alg = HSLDDRK64()
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
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
# test whether aliasing u0 is bad
new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
@test sol_old[end] ≈ sol_new[end]


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
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
# test whether aliasing u0 is bad
new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
@test sol_old[end] ≈ sol_new[end]


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
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
# test whether aliasing u0 is bad
new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
@test sol_old[end] ≈ sol_new[end]


alg = DGLDDRK84_F()
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
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
# test whether aliasing u0 is bad
new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
@test sol_old[end] ≈ sol_new[end]


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
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
# test whether aliasing u0 is bad
new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
@test sol_old[end] ≈ sol_new[end]


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
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
integ = init(prob_ode_large, alg, dt=1.e-2, save_start=false, save_end=false, save_everystep=false, alias_u0=true)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3
# test whether aliasing u0 is bad
new_prob_ode_nonlinear_inplace = ODEProblem(prob_ode_nonlinear_inplace.f,[1.],(0.,0.5))
sol_old = solve(prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false)
sol_new = solve(new_prob_ode_nonlinear_inplace, alg, dt=1.e-4, save_everystep=false, save_start=false, alias_u0=true)
@test sol_old[end] ≈ sol_new[end]


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

println("Methods from Carpenter, Kennedy, Lewis (2000)")

alg = CKLLSRK43_2()
dts = 1 ./ 2 .^(8:-1:4)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+1 atol=testTol    # This scheme has kinear order of 4
end
for prob in test_problems_nonlinear
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


alg = CKLLSRK54_3C()
dts = 1 ./ 2 .^(8:-1:4)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ 1 atol=testTol          # The CI plot is linear but the evaluated order is 1
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
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


alg = CKLLSRK95_4S()
dts = 1 ./ 2 .^(6:-1:2)
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


alg = CKLLSRK95_4C()
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


alg = CKLLSRK95_4M()
dts = 1 ./ 2 .^(6:-1:2)
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


alg = CKLLSRK54_3C_3R()
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


alg = CKLLSRK54_3M_3R()
dts = 1 ./ 2 .^(7:-1:3)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg)+1 atol=testTol
end
for prob in test_problems_nonlinear
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


alg = CKLLSRK54_3N_3R()
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


alg = CKLLSRK85_4C_3R()
dts = 1 ./ 2 .^(6:-1:2)
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


alg = CKLLSRK85_4M_3R()
dts = 1.6 ./ 2 .^(6:-1:2)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
dts = 1 ./ 2 .^(6:-1:2)
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
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


alg = CKLLSRK85_4P_3R()
dts = 1 ./ 2 .^(6:-1:2)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
dts = 1.0 ./ 10.0 .^(1.3:-.22:0.3)
for prob in test_problems_linear
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


alg = CKLLSRK54_3N_4R()
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


alg = CKLLSRK54_3M_4R()
dts = 1 ./ 1.9 .^(7:-1:3)
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
dts = 1 ./ 1.8 .^(8:-1:4)
for prob in test_problems_linear
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


alg = CKLLSRK65_4M_4R()
dts = 1 ./ 2 .^(6:-1:2)
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


alg = CKLLSRK85_4FM_4R()
dts = (1) ./ 2 .^(10:-1:6)
for prob in test_problems_only_time
	sim = test_convergence(dts, prob, alg)
	@test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) + 1 atol=testTol
end
for prob in test_problems_nonlinear
	sim = test_convergence(dts, prob, alg)
	@test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
end
for prob in test_problems_linear
	sim = test_convergence(dts, prob, alg)
	@test_broken sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
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


alg = CKLLSRK75_4M_5R()
dts = (1) ./ 1.9 .^(6:-1:2)
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


println("Methods from Parsani, Ketcheson, Deconinck (2013)")

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
dts = 1 ./ 2 .^(10:-1:8)
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
