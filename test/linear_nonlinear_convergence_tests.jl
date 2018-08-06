using OrdinaryDiffEq, Test, DiffEqDevTools, DiffEqOperators, Random
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
  for Alg in [GenericIIF1,GenericIIF2,LawsonEuler,NorsettEuler,ETDRK2,ETDRK3,ETDRK4,HochOst4,Exprb32,Exprb43,ETD2]
    sim  = test_convergence(dts,prob,Alg())
    if Alg in [Exprb32, Exprb43]
      @test_broken abs(sim.𝒪est[:l2] - alg_order(Alg())) < 0.2
    else
      @test abs(sim.𝒪est[:l2] - alg_order(Alg())) < 0.2
    end
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

  dts = 1 ./2 .^(8:-1:4) #14->7 good plot
  for Alg in [GenericIIF1,GenericIIF2,LawsonEuler,NorsettEuler,ETDRK2,ETDRK3,ETDRK4,HochOst4,Exprb32,Exprb43,ETD2]
    sim  = test_convergence(dts,prob,Alg())
    if Alg in [Exprb32, Exprb43]
      @test_broken abs(sim.𝒪est[:l2] - alg_order(Alg())) < 0.1
    else
      @test abs(sim.𝒪est[:l2] - alg_order(Alg())) < 0.1
    end
  end
  sim  = test_convergence(dts,prob,ETDRK4(),dense_errors=true)
  @test abs(sim.𝒪est[:l2]-4) < 0.1
  @test abs(sim.𝒪est[:L2]-4) < 0.1
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
      @test_broken abs(sim.𝒪est[:l2] - alg_order(Alg())) < 0.1
    else
      @test abs(sim.𝒪est[:l2] - alg_order(Alg())) < 0.1
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
  srand(0); u0 = rand(2); tspan = (0.0, 1.0)
  prob = ODEProblem(fun, u0, tspan)
  # Setup approximate solution
  test_setup = Dict(:alg=>Vern9(), :reltol=>1e-16, :abstol=>1e-16)
  # Convergence simulation
  dts = 1 ./2 .^(7:-1:4)
  Algs = [Exp4, EPIRK4s3A, EPIRK4s3B, EPIRK5s3, EXPRB53s3, EPIRK5P1, EPIRK5P2]
  for Alg in Algs
    sim = analyticless_test_convergence(dts, prob, Alg(adaptive_krylov=false), test_setup)
    if Alg == EPIRK5s3
      @test_broken abs(sim.𝒪est[:l2] - alg_order(Alg())) < 0.1
    else
      @test abs(sim.𝒪est[:l2] - alg_order(Alg())) < 0.1
    end
  end
end
