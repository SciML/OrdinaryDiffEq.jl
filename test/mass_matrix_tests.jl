using OrdinaryDiffEq, Test

# TODO: clean up
@testset "Mass Matrix Accuracy Tests" begin
  mm_A = [-2.0 1 4
           4 -2 1
           2 1 3]
  mm_b = mm_A*ones(3)
  function mm_f(du,u,p,t)
    A_mul_B!(du,mm_A,u)
    tmp = t*mm_b
    du .+= tmp
  end
  mm_f(::Type{Val{:analytic}},u0,p,t) = @. 2ones(3)*exp(t) - t - 1
  mm_g(du,u,p,t) = du .= u + t
  mm_g(::Type{Val{:analytic}},u0,p,t) = @. 2ones(3)*exp(t) - t - 1
  prob2 = ODEProblem(mm_g,ones(3),(0.0,1.0))
  prob = ODEProblem(mm_f,ones(3),(0.0,1.0),mass_matrix=mm_A)

  ######################################### Test each method for exactness

  sol = solve(prob  ,Rosenbrock23())
  sol2 = solve(prob2,Rosenbrock23())

  @test norm(sol .- sol2) ≈ 0 atol=1e-11

  sol = solve(prob, Rosenbrock32())
  sol2 = solve(prob2,Rosenbrock32())

  @test norm(sol .- sol2) ≈ 0 atol=1e-11

  sol = solve(prob,  ROS3P())
  sol2 = solve(prob2,ROS3P())

  @test norm(sol .- sol2) ≈ 0 atol=1e-11

  sol = solve(prob,  Rodas3())
  sol2 = solve(prob2,Rodas3())

  @test norm(sol .- sol2) ≈ 0 atol=1e-11

  sol = solve(prob,  RosShamp4())
  sol2 = solve(prob2,RosShamp4())

  @test norm(sol .- sol2) ≈ 0 atol=1e-10

  sol = solve(prob,  Rodas4())
  sol2 = solve(prob2,Rodas4())

  @test norm(sol .- sol2) ≈ 0 atol=1e-9

  sol = solve(prob,  Rodas5())
  sol2 = solve(prob2,Rodas5())

  @test norm(sol .- sol2) ≈ 0 atol=1e-7

  sol = solve(prob,  ImplicitEuler(),dt=1/10,adaptive=false)
  sol2 = solve(prob2,ImplicitEuler(),dt=1/10,adaptive=false)

  @test norm(sol .- sol2) ≈ 0 atol=1e-7

  sol = solve(prob,  ImplicitMidpoint(),dt=1/10)
  sol2 = solve(prob2,ImplicitMidpoint(),dt=1/10)

  @test norm(sol .- sol2) ≈ 0 atol=1e-7
end

#=

sol = solve(prob,  Trapezoid())
sol2 = solve(prob2,Trapezoid())

@test norm(sol .- sol2) ≈ 0 atol=1e-7

sol = solve(prob,  TRBDF2())
sol2 = solve(prob2,TRBDF2())

@test norm(sol .- sol2) ≈ 0 atol=1e-9

#sol = solve(prob, SDIRK2())
#sol2 = solve(prob2, SDIRK2())

#@test norm(sol .- sol2) ≈ 0 atol=1e-11

sol = solve(prob,   TRBDF2(),adaptive=false,dt=1/10)
sol2 = solve(prob2, TRBDF2(),adaptive=false,dt=1/10)

=#

# Singular mass matrices

@testset "Mass Matrix Tests with Singular Mass Matrices" begin
  function f!(du, u, p, t)
      du[1] = u[2]
      du[2] = u[2] - 1.
      return
  end

  u0 = [0.,1.]
  tspan = (0.0, 1.0)

  M = zeros(2,2)
  M[1,1] = 1.

  m_ode_prob = ODEProblem(f!, u0, tspan, mass_matrix=M)
  @test_nowarn sol = solve(m_ode_prob, Rosenbrock23())

  M = [0.637947  0.637947
       0.637947  0.637947]

  inv(M) # not caught as singular

  function f2!(du, u, p, t)
      du[1] = u[2]
      du[2] = u[1]
      return
  end
  u0 = zeros(2)

  m_ode_prob = ODEProblem(f2!, u0, tspan, mass_matrix=M)
  @test_nowarn sol = solve(m_ode_prob, Rosenbrock23())
end
