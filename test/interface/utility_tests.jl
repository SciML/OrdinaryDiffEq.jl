using OrdinaryDiffEq: WOperator, calc_W, calc_W!
using OrdinaryDiffEq, LinearAlgebra, SparseArrays, Random, Test, DiffEqOperators

@testset "WOperator" begin
  Random.seed!(123)
  y = zeros(2); b = rand(2)
  mm = rand(2, 2)
  for _J in [rand(2, 2), Diagonal(rand(2))]
    _Ws = [-mm + 2.0 * _J, -mm/2.0 + _J]
    for inplace in (true, false), (_W, W_transform) in zip(_Ws, [false, true])
      W = WOperator{inplace}(mm, 2.0, DiffEqArrayOperator(_J), b, transform=W_transform)
      @test convert(AbstractMatrix,W) ≈ _W
      @test W * b ≈ _W * b
      mul!(y, W, b); @test y ≈ _W * b
    end
  end
end

@testset "calc_W and calc_W!" begin
  A = [-1.0 0.0; 0.0 -0.5]; mm = [2.0 0.0; 0.0 1.0]
  u0 = [1.0, 1.0]; tmp = zeros(2)
  tspan = (0.0,1.0); dt = 0.01; dtgamma = 0.5dt
  concrete_W = -mm + dtgamma * A

  # Out-of-place
  fun = ODEFunction((u,p,t) -> A*u;
                    mass_matrix=mm,
                    jac=(u,p,t) -> A)
  integrator = init(ODEProblem(fun,u0,tspan), ImplicitEuler(); adaptive=false, dt=dt)
  W = calc_W(integrator, integrator.cache.nlsolver.cache, dtgamma, false)
  @test convert(AbstractMatrix, W) == concrete_W
  @test W \ u0 ≈ concrete_W \ u0

  # In-place
  fun = ODEFunction((du,u,p,t) -> mul!(du,A,u);
                    mass_matrix=mm,
                    jac_prototype=DiffEqArrayOperator(A))
  integrator = init(ODEProblem(fun,u0,tspan), ImplicitEuler(); adaptive=false, dt=dt)
  calc_W!(integrator.cache.nlsolver.cache.W, integrator, integrator.cache.nlsolver, integrator.cache, dtgamma, false)
  @test convert(AbstractMatrix, integrator.cache.nlsolver.cache.W) == concrete_W
  ldiv!(tmp, lu!(integrator.cache.nlsolver.cache.W), u0); @test tmp == concrete_W \ u0
end

@testset "Implicit solver with lazy W" begin
  A = sparse([-1.0 0.0; 0.0 -0.5])
  mm = sparse([2.0 0.0; 0.0 1.0])
  u0 = [1.0, 1.0]; tspan = (0.0,1.0)

  _f = (u,p,t) -> t*(A*u); _f_ip = (du,u,p,t) -> lmul!(t,mul!(du,A,u))
  fun1 = ODEFunction(_f; mass_matrix=mm)
  fun2 = ODEFunction(_f; mass_matrix=mm, jac=(u,p,t) -> t*A)
  fun1_ip = ODEFunction(_f_ip; mass_matrix=mm)
  fun2_ip = ODEFunction(_f_ip; mass_matrix=mm,
  jac_prototype=DiffEqArrayOperator(similar(A); update_func=(J,u,p,t) -> (J .= t .* A; J)))

  for Alg in [ImplicitEuler, Rosenbrock23, Rodas5]
    println(Alg)
    sol1 = solve(ODEProblem(fun1,u0,tspan), Alg(); adaptive=false, dt=0.01)
    sol2 = solve(ODEProblem(fun2,u0,tspan), Alg(); adaptive=false, dt=0.01)
    @test sol1(1.0) ≈ sol2(1.0)

    sol1_ip = solve(ODEProblem(fun1_ip,u0,tspan), Alg(); adaptive=false, dt=0.01)
    sol2_ip = solve(ODEProblem(fun2_ip,u0,tspan), Alg(linsolve=LinSolveFactorize(lu)); adaptive=false, dt=0.01)
    @test sol1_ip(1.0) ≈ sol2_ip(1.0) atol=1e-5
  end
end
