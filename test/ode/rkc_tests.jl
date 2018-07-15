using OrdinaryDiffEq, DiffEqDevTools, Test, LinearAlgebra
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
using OrdinaryDiffEq: maxeig!
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear
probArr = Vector{ODEProblem}(undef, 2)
probArr[1] = prob_ode_linear
probArr[2] = prob_ode_2Dlinear

srand(123)
@testset "Power Iteration of Runge-Kutta-Chebyshev Tests" begin
  for i in 1:10
    A = rand(200,200)
    test_f(u,p,t) = A*u
    prob = ODEProblem(test_f, rand(200), (0,1.))
    integrator = init(prob, ROCK2())
    eigm = maximum(abs.(eigvals(A)))
    maxeig!(integrator, integrator.cache)
    eigest = integrator.eigen_est
    @test eigest ≈ eigm atol=0.22eigm
  end
end

@testset "Runge-Kutta-Chebyshev Convergence Tests" begin
  for prob in probArr
    sim2 = test_convergence(dts,prob,ROCK2())
    @test abs(sim2.𝒪est[:l∞]-2) < testTol
  end
end
