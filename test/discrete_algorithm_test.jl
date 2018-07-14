using OrdinaryDiffEq, Test, Statistics
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_2Dlinear, prob_ode_linear

@testset "Scalar Discrete Problem" begin
  prob = DiscreteProblem(0.5,(0.0,1.0))
  sol =solve(prob,FunctionMap())
  @test sol[1] == sol[end]
  @test sol(0.5:0.1:0.7) == [sol[1],sol[1],sol[1]]

  sol =solve(prob_ode_linear,FunctionMap())
  @test sol[end] ≈ .505
  sol =solve(prob_ode_linear,FunctionMap(scale_by_time=true),dt=1/4)
  sol2 =solve(prob_ode_linear,Euler(),dt=1/4)
  @test sol[end] == sol2[end]
  @test sol(0.53) != sol2(0.53)
end

@testset "Array Discrete Problem" begin
  prob2 = DiscreteProblem(rand(4,2),(0.0,1.0))
  sol =solve(prob2,FunctionMap())
  @test sol[1] == sol[end]
  @test sol(0.5) == sol[1]

  @test_nowarn sol =solve(prob_ode_2Dlinear,FunctionMap())
  @test_nowarn sol2 =solve(prob_ode_2Dlinear,Euler(),dt=1)
  sol =solve(prob_ode_2Dlinear,FunctionMap(scale_by_time=true),dt=1/4)
  sol2 =solve(prob_ode_2Dlinear,Euler(),dt=1/4)
  @test sol[end] == sol2[end]
  @test sol(0.35) != sol2(0.53)
end

@testset "Discrete Problem (Henon Map)" begin
  function henon_map!(u_next, u, _p, t)
    u_next[1] = 1 + u[2] - 1.4 * u[1]^2
    u_next[2] = 0.3 * u[1]
  end

  prob = DiscreteProblem(henon_map!, [0.5, 0.5], (0, 100)) # Integer time
  sol = solve(prob, FunctionMap())
  @test eltype(sol.t) <: Integer
end
