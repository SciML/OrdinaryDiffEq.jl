using OrdinaryDiffEq, DiffEqDevTools, Test, LinearAlgebra, Random
using OrdinaryDiffEq: maxeig!
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear
probArr = Vector{ODEProblem}(undef, 2)
probArr[1] = prob_ode_linear
probArr[2] = prob_ode_2Dlinear

@testset "Power Iteration of Runge-Kutta-Chebyshev Tests" begin
  Random.seed!(123)
  eigen_est = (integrator) -> integrator.eigen_est = 1.5e2
  for iip in [true, false], alg in [ROCK4(), ROCK4(eigen_est=eigen_est)]
    println(typeof(alg))
    A = randn(20,20)
    test_f(u,p,t) = A*u
    test_f(du,u,p,t) = mul!(du, A, u)
    prob = ODEProblem{iip}(test_f, randn(20), (0,1.))
    integrator = init(prob, alg)
    eigm = maximum(abs.(eigvals(A)))
    maxeig!(integrator, integrator.cache)
    eigest = integrator.eigen_est
    @test eigest ≈ eigm rtol=0.1eigm

    A = A - 1e2I
    test_stiff(u,p,t) = A*u
    test_stiff(du,u,p,t) = mul!(du, A, u)
    prob = ODEProblem{iip}(test_stiff, ones(20), (0,1.))
    @test_nowarn solve(prob, alg)
  end

  Random.seed!(123)
  for iip in [true, false], Alg in [IRKC]
    alg = Alg()
    println(typeof(alg))
    A = randn(20,20)
    B = randn(20,20)
    test_f1 = !iip ? (u,p,t) -> A*u : (du,u,p,t) -> mul!(du, A, u)
    test_f2 = !iip ? (u,p,t) -> B*u : (du,u,p,t) -> mul!(du, B, u)
    ff_split = SplitFunction{iip}(test_f1,test_f2)
    prob = SplitODEProblem{iip}(ff_split, randn(20,1), (0.0,1.))
    integrator = init(prob, alg)
    eigm = maximum(abs.(eigvals(A)))
    maxeig!(integrator, integrator.cache)
    eigest = integrator.eigen_est
    @test eigest ≈ eigm rtol=0.1eigm

    A = A - 1e2I
    test_f1 = !iip ? (u,p,t) -> A*u : (du,u,p,t) -> mul!(du, A, u)
    prob = SplitODEProblem{iip}(SplitFunction{iip}(test_f1,test_f2), ones(20), (0.0, 1.0))
    @test_nowarn solve(prob, alg)
  end
end

@testset "Runge-Kutta-Chebyshev Convergence Tests" begin
  dts = 1 .//2 .^(8:-1:4)
  testTol = 0.1
  for prob in probArr
    println("ROCK2")
    #default ROCK2
    sim = test_convergence(dts,prob,ROCK2())
    @test sim.𝒪est[:l∞] ≈ 2 atol=testTol
    #testing ROCK2 for different minimum stages to insure that the constatns are right
    sim = test_convergence(dts,prob,ROCK2(min_stages=5))
    @test sim.𝒪est[:l∞] ≈ 2 atol=testTol
    sim = test_convergence(dts,prob,ROCK2(min_stages=10))
    @test sim.𝒪est[:l∞] ≈ 2 atol=testTol
    sim = test_convergence(dts,prob,ROCK2(min_stages=21))
    @test sim.𝒪est[:l∞] ≈ 2 atol=testTol
    #default ROCK4
    println("ROCK4")
    sim = test_convergence(dts,prob,ROCK4())
    @test sim.𝒪est[:l∞] ≈ 4 atol=testTol
    #testing ROCK4 for different minimum stages to insure that the constatns are right
    sim = test_convergence(dts,prob,ROCK4(min_stages=6))
    @test sim.𝒪est[:l∞] ≈ 4 atol=testTol
    sim = test_convergence(dts,prob,ROCK4(min_stages=10))
    @test sim.𝒪est[:l∞] ≈ 4 atol=testTol
    sim = test_convergence(dts,prob,ROCK4(min_stages=21))
    @test sim.𝒪est[:l∞] ≈ 4 atol=testTol

    println("ROCKC")
    sim = test_convergence(dts,prob,RKC())
    @test sim.𝒪est[:l∞] ≈ 2 atol=testTol
    println("SERK2")
    sim = test_convergence(dts,prob,SERK2())
    @test sim.𝒪est[:l∞] ≈ 2 atol=testTol
    println("ESERK4")
    sim = test_convergence(dts,prob,ESERK4())
    @test sim.𝒪est[:l∞] ≈ 4 atol=testTol
  end
  dts = 1 .//2 .^(6:-1:2)
  for prob in probArr
    println("ESERK5")
    sim = test_convergence(dts,prob,ESERK5())
    @test sim.𝒪est[:l∞] ≈ 5 atol=testTol
  end
end
