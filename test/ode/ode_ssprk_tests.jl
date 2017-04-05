using OrdinaryDiffEq, DiffEqProblemLibrary, DiffEqDevTools, Base.Test

dts = 1.//2.^(8:-1:4)
testTol = 0.2

prob_ode_sin = ODETestProblem((t,u)->cos(t), 0., (t,u0)->sin(t))
prob_ode_sin_inplace = ODETestProblem((t,u,du)->du[1]=cos(t), [0.], (t,u0)->[sin(t)])

test_problems = [prob_ode_sin, prob_ode_sin_inplace,
                 prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear]

for prob in test_problems
  for alg in [SSPRK22()]
    sim = test_convergence(dts, prob, alg)
    @test abs(sim.𝒪est[:final]-OrdinaryDiffEq.alg_order(alg)) < testTol
  end
end
