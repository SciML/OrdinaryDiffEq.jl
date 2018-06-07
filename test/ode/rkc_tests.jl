using OrdinaryDiffEq, DiffEqDevTools, DiffEqProblemLibrary, Base.Test
using OrdinaryDiffEq: maxeig!

srand(123)
@testset "Power Iteration of Runge-Kutta-Chebyshev" begin
  for i in 1:10
    A = rand(200,200)
    f(u,p,t) = A*u
    prob = ODEProblem(f, rand(200), (0,1.))
    integrator = init(prob, ROCK2())
    eigm = maximum(abs.(eigvals(A)))
    maxeig!(integrator, integrator.cache)
    eigest = integrator.eigen_est
    @test eigest ≈ eigm atol=0.22eigm
  end
end
