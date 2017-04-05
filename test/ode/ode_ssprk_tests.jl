using DiffEqBase, OrdinaryDiffEq, DiffEqProblemLibrary, DiffEqDevTools, Base.Test

dts = 1.//2.^(8:-1:4)
testTol = 0.25

prob_ode_sin = ODETestProblem((t,u)->cos(t), 0., (t,u0)->sin(t))
prob_ode_sin_inplace = ODETestProblem((t,u,du)->du[1]=cos(t), [0.], (t,u0)->[sin(t)])

prob_ode_nonlinear = ODETestProblem( (t,u)->sin(u), 1.,
                                     (t,u0)->2*acot(exp(-t)*cot(0.5)), (0.,0.5) )
prob_ode_nonlinear_inplace = ODETestProblem( (t,u,du)->du[1]=sin(u[1]), [1.],
                                             (t,u0)->[2*acot(exp(-t)*cot(0.5))], (0.,0.5) )

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
