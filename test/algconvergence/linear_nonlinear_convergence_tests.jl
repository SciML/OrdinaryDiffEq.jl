using OrdinaryDiffEq, Test, DiffEqDevTools, DiffEqOperators, Random, LinearAlgebra
using OrdinaryDiffEq: alg_order

@testset "Caching Out-of-place" begin
  println("Caching Out-of-place")
  μ = 1.01
  linnonlin_f2 = (u,p,t) -> μ * u
  linnonlin_f1 = DiffEqScalar(μ)
  linnonlin_fun = SplitFunction(linnonlin_f1, linnonlin_f2; analytic=(u0,p,t)->u0.*exp.(2μ*t))
  prob = SplitODEProblem(linnonlin_fun,1/2,(0.0,1.0))

  Random.seed!(100)
  dts = 1 ./2 .^(7:-1:4) #14->7 good plot
  for Alg in [LawsonEuler,NorsettEuler,ETDRK2,ETDRK3,ETDRK4,HochOst4,ETD2,KenCarp3,CFNLIRK3,CFNLIRK4]
    sim  = test_convergence(dts,prob,Alg())
    @test sim.𝒪est[:l2] ≈ alg_order(Alg()) atol=0.2
  end
  # Dense test
  sim  = test_convergence(dts,prob,ETDRK4(),dense_errors=true)
  @test abs(sim.𝒪est[:l2]-4) < 0.2
end

@testset "Caching Inplace" begin
  println("Caching Inplace")
  μ = 1.01
  u0 = rand(2)
  A = [2.0 -1.0; -1.0 2.0]
  linnonlin_f1 = DiffEqArrayOperator(A)
  linnonlin_f2 = (du,u,p,t) -> du .= μ .* u
  linnonlin_fun_iip = SplitFunction(linnonlin_f1,linnonlin_f2;analytic=(u0,p,t)->exp((A+μ*I)*t)*u0)
  prob = SplitODEProblem(linnonlin_fun_iip,u0,(0.0,1.0))

  dts = 1 ./2 .^(7:-1:4) #14->7 good plot
  for Alg in [LawsonEuler(),NorsettEuler(),ETDRK2(),ETDRK3(),ETDRK4(),HochOst4(),ETD2(),KenCarp3(linsolve=LinSolveGMRES(tol=1e-6)),CFNLIRK3(),CFNLIRK4()]
    sim  = test_convergence(dts,prob,Alg)
    @test sim.𝒪est[:l2] ≈ alg_order(Alg) atol=0.17
  end
  sim  = test_convergence(dts,prob,ETDRK4(),dense_errors=true)
  @test sim.𝒪est[:l2] ≈  4 atol=0.1
  @test sim.𝒪est[:L2] ≈ 4 atol=0.1
end

@testset "EPIRK Out-of-place" begin
  println("EPIRK Out-of-place")
  # Setup nonlinear problem
  A = [-2.0 1.0; 1.0 -2.0]
  f = (u,p,t) -> A*u - u.^3
  jac = (u,p,t) -> A - [3u[1]^2 0.0; 0.0 3u[2]^2]
  fun = ODEFunction(f; jac=jac)
  Random.seed!(0); u0 = rand(2); tspan = (0.0, 1.0)
  prob = ODEProblem(fun, u0, tspan)
  # Setup approximate solution
  test_setup = Dict(:alg=>Vern9(), :reltol=>1e-16, :abstol=>1e-16)
  # Convergence simulation
  dts = 1 ./2 .^(7:-1:4)
  Algs = [Exp4, EPIRK4s3A, EPIRK4s3B, EPIRK5s3, EXPRB53s3, EPIRK5P1, EPIRK5P2]
  for Alg in Algs
    sim = analyticless_test_convergence(dts, prob, Alg(adaptive_krylov=false), test_setup)
    if Alg == EPIRK5s3
      @test_broken sim.𝒪est[:l2] ≈ alg_order(Alg()) atol=0.1
    else
      @test sim.𝒪est[:l2] ≈ alg_order(Alg()) atol=0.1
    end
  end
end

@testset "EPIRK Inplace" begin
  println("EPIRK Inplace")
  # Setup nonlinear problem
  A = [-2.0 1.0; 1.0 -2.0]
  f = (du,u,p,t) -> (mul!(du, A, u); du .-= u.^3)
  jac_update = (J,u,p,t) -> (copyto!(J,A); J[1,1] -= 3u[1]^2; J[2,2] -= 3u[2]^2)
  jac_prototype = DiffEqArrayOperator(zeros(2,2); update_func=jac_update)
  fun = ODEFunction(f; jac_prototype=jac_prototype)
  Random.seed!(0); u0 = rand(2); tspan = (0.0, 1.0)
  prob = ODEProblem(fun, u0, tspan)
  # Setup approximate solution
  test_setup = Dict(:alg=>Vern9(), :reltol=>1e-16, :abstol=>1e-16)
  # Convergence simulation
  dts = 1 ./2 .^(7:-1:4)
  Algs = [Exp4, EPIRK4s3A, EPIRK4s3B, EPIRK5s3, EXPRB53s3, EPIRK5P1, EPIRK5P2]
  for Alg in Algs
    sim = analyticless_test_convergence(dts, prob, Alg(adaptive_krylov=false), test_setup)
    if Alg == EPIRK5s3
      @test_broken sim.𝒪est[:l2] ≈ alg_order(Alg()) atol=0.1
    else
      @test sim.𝒪est[:l2] ≈ alg_order(Alg()) atol=0.1
    end
  end
end

@testset "ExpRK Autodiff" begin
  println("ExpRK Autodiff")
  # Setup nonlinear problem (without explicit jacobian)
  A = [-2.0 1.0; 1.0 -2.0]
  f = (u,p,t) -> A*u - u.^3
  f_ip = (du,u,p,t) -> (mul!(du, A, u); du .-= u.^3)
  fun = ODEFunction(f)
  fun_ip = ODEFunction(f_ip)
  Random.seed!(0); u0 = rand(2); tspan = (0.0, 1.0)
  prob = ODEProblem(fun, u0, tspan)
  prob_ip = ODEProblem(fun_ip, u0, tspan)
  # Setup approximate solution
  test_setup = Dict(:alg=>Vern9(), :reltol=>1e-16, :abstol=>1e-16)
  # Convergence simulation
  dts = 1 ./2 .^(7:-1:4)
  sim = analyticless_test_convergence(dts, prob, HochOst4(krylov=true), test_setup)
  @test sim.𝒪est[:l2] ≈ 4 atol=0.1
  sim = analyticless_test_convergence(dts, prob_ip, HochOst4(krylov=true), test_setup)
  @test sim.𝒪est[:l2] ≈ 4 atol=0.1
  sim = analyticless_test_convergence(dts, prob, EPIRK5P1(adaptive_krylov=false), test_setup)
  @test sim.𝒪est[:l2] ≈ 5 atol=0.1
  sim = analyticless_test_convergence(dts, prob_ip, EPIRK5P1(adaptive_krylov=false), test_setup)
  @test sim.𝒪est[:l2] ≈ 5 atol=0.1
end

@testset "Adaptive Exprb Out-of-place" begin
  println("Adaptive Exprb Out-of-place")
  # Setup nonlinear problem
  A = [-2.0 1.0; 1.0 -2.0]
  f = (u,p,t) -> A*u - u.^3
  jac = (u,p,t) -> A - [3u[1]^2 0.0; 0.0 3u[2]^2]
  fun = ODEFunction(f; jac=jac)
  Random.seed!(0); u0 = rand(2); tspan = (0.0, 1.0)
  prob = ODEProblem(fun, u0, tspan)
  # Setup approximate solution
  test_setup = Dict(:alg=>Vern9(), :reltol=>1e-16, :abstol=>1e-16)
  # Convergence simulation
  dts = 1 ./2 .^(7:-1:4)
  Algs = [Exprb32, Exprb43]
  for Alg in Algs
    sim = analyticless_test_convergence(dts, prob, Alg(), test_setup)
    @test sim.𝒪est[:l2] ≈ alg_order(Alg()) atol=0.1
  end
end

@testset "Adaptive Exprb Inplace" begin
  println("Adaptive Exprb Inplace")
  # Setup nonlinear problem
  A = [-2.0 1.0; 1.0 -2.0]
  f = (du,u,p,t) -> (mul!(du, A, u); du .-= u.^3)
  jac_update = (J,u,p,t) -> (copyto!(J,A); J[1,1] -= 3u[1]^2; J[2,2] -= 3u[2]^2)
  jac_prototype = DiffEqArrayOperator(zeros(2,2); update_func=jac_update)
  fun = ODEFunction(f; jac_prototype=jac_prototype)
  Random.seed!(0); u0 = rand(2); tspan = (0.0, 1.0)
  prob = ODEProblem(fun, u0, tspan)
  # Setup approximate solution
  test_setup = Dict(:alg=>Vern9(), :reltol=>1e-16, :abstol=>1e-16)
  # Convergence simulation
  dts = 1 ./2 .^(7:-1:4)
  Algs = [Exprb32, Exprb43]
  for Alg in Algs
    sim = analyticless_test_convergence(dts, prob, Alg(), test_setup)
    @test sim.𝒪est[:l2] ≈ alg_order(Alg()) atol=0.1
  end
end
