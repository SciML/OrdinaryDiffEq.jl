using OrdinaryDiffEq
using Test

let prob = ODEProblem((du,u,p,t)->du.=u,ones(2),(0, 1.))
  # issue 635
  for alg in (ImplicitEuler(), Trapezoid())
    @test_logs (:warn, "Adaptive $(nameof(typeof(alg))) method is not safe to use with save_everystep=false due to its error estimator.") solve(prob, alg, save_everystep=false)
    @test_nowarn solve(prob, alg, adaptive=false, dt=0.1, save_everystep=false)
    @test_nowarn solve(prob, alg)
  end

  # saveat
  @test_logs (:warn, "Dense output is incompatible with saveat. Please use the SavingCallback from the Callback Library to mix the two behaviors.") solve(prob, Tsit5(), dense=true, saveat=0.1)

  @test_throws ErrorException solve(prob, Euler())
end

# mass matrix
@test_throws ErrorException solve(ODEProblem(ODEFunction((du,u,p,t)->du.=u, mass_matrix=rand(2, 2)),ones(2),(0, 1.)), Tsit5())
