@testset "Runge-Kutta-Chebyshev Tests" begin
using OrdinaryDiffEq, DiffEqDevTools, Test, LinearAlgebra, Random
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
using OrdinaryDiffEq: maxeig!
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear
probArr = Vector{ODEProblem}(undef, 2)
probArr[1] = prob_ode_linear
probArr[2] = prob_ode_2Dlinear

Random.seed!(123)
@testset "Power Iteration of Runge-Kutta-Chebyshev Tests" begin
  for i in 1:10, iip in [true, false]
    A = randn(20,20)
    test_f(u,p,t) = A*u
    test_f(du,u,p,t) = mul!(du, A, u)
    prob = ODEProblem{iip}(test_f, randn(20), (0,1.))
    integrator = init(prob, ROCK2())
    eigm = maximum(abs.(eigvals(A)))
    maxeig!(integrator, integrator.cache)
    eigest = integrator.eigen_est
    @test eigest ≈ eigm rtol=0.1eigm
  end
end

@testset "Runge-Kutta-Chebyshev Convergence Tests" begin
  dts = 1 .//2 .^(8:-1:4)
  testTol = 0.1
  for prob in probArr
    sim2 = test_convergence(dts,prob,ROCK2())
    @test abs(sim2.𝒪est[:l∞]-2) < testTol
  end
end
end
