using OrdinaryDiffEq, DiffEqProblemLibrary, DiffEqDevTools, Base.Test

dts = 1.//2.^(8:-1:4)
testTol = 0.2

for prob in [prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear]
  for alg in [SSPRK22()]
    sim = test_convergence(dts, prob, alg)
    @test abs(sim.𝒪est[:final]-OrdinaryDiffEq.alg_order(alg)) < testTol
  end
end
