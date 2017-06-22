using DiffEqBase, OrdinaryDiffEq, DiffEqProblemLibrary, DiffEqDevTools, Base.Test

dts = 1.//2.^(8:-1:4)
testTol = 0.25

f = (t,u)->cos(t)
(p::typeof(f))(::Type{Val{:analytic}},t,u0) = sin(t)
prob_ode_sin = ODEProblem(f, 0.,(0.0,1.0))

f = (t,u,du)->du[1]=cos(t)
(p::typeof(f))(::Type{Val{:analytic}},t,u0) = [sin(t)]
prob_ode_sin_inplace = ODEProblem(f, [0.], (0.0,1.0))

f = (t,u)->sin(u)
(p::typeof(f))(::Type{Val{:analytic}},t,u0) = 2*acot(exp(-t)*cot(0.5))
prob_ode_nonlinear = ODEProblem(f, 1.,(0.,0.5))

f = (t,u,du)->du[1]=sin(u[1])
(p::typeof(f))(::Type{Val{:analytic}},t,u0) = [2*acot(exp(-t)*cot(0.5))]
prob_ode_nonlinear_inplace = ODEProblem(f,[1.],(0.,0.5))

const linear_bigα = parse(BigFloat,"1.01")
f_2dlinearbig = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = linear_bigα*u[i]
  end
end
(f::typeof(f_2dlinearbig))(::Type{Val{:analytic}},t,u0) = u0*exp.(1.01*t)
prob_ode_bigfloat2Dlinear = ODEProblem(f_2dlinearbig,map(BigFloat,rand(4,2)).*ones(4,2)/2,(0.0,1.0))


test_problems_only_time = [prob_ode_sin, prob_ode_sin_inplace]
test_problems_linear = [prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear]
test_problems_nonlinear = [prob_ode_nonlinear, prob_ode_nonlinear_inplace]


alg = SSPRK22()
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test abs(sim.𝒪est[:final]-OrdinaryDiffEq.alg_order(alg)) < testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test abs(sim.𝒪est[:final]-OrdinaryDiffEq.alg_order(alg)) < testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test abs(sim.𝒪est[:final]-OrdinaryDiffEq.alg_order(alg)) < testTol
end


alg = SSPRK33()
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  # This corresponds to Simpson's rule; due to symmetric quadrature nodes,
  # it is of degree 4 instead of 3, as would be expected.
  @test abs(sim.𝒪est[:final]-1-OrdinaryDiffEq.alg_order(alg)) < testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test abs(sim.𝒪est[:final]-OrdinaryDiffEq.alg_order(alg)) < testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test abs(sim.𝒪est[:final]-OrdinaryDiffEq.alg_order(alg)) < testTol
end


alg = SSPRK104()
for prob in test_problems_only_time
  sim = test_convergence(dts, prob, alg)
  @test abs(sim.𝒪est[:final]-OrdinaryDiffEq.alg_order(alg)) < testTol
end
for prob in test_problems_linear
  sim = test_convergence(dts, prob, alg)
  @test abs(sim.𝒪est[:final]-OrdinaryDiffEq.alg_order(alg)) < testTol
end
for prob in test_problems_nonlinear
  sim = test_convergence(dts, prob, alg)
  @test abs(sim.𝒪est[:final]-OrdinaryDiffEq.alg_order(alg)) < testTol
end
