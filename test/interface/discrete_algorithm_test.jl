using OrdinaryDiffEq, Test, Statistics
import ODEProblemLibrary: prob_ode_2Dlinear, prob_ode_linear

@testset "Scalar Discrete Problem" begin
    prob = DiscreteProblem(0.5, (0.0, 1.0))
    sol = solve(prob, FunctionMap())
    @test sol[1] == sol[end]
    @test sol(0.5:0.1:0.7) == [sol[1], sol[1], sol[1]]

    sol = solve(prob_ode_linear, FunctionMap())
    @test sol[end] ≈ 0.505
    sol = solve(prob_ode_linear, FunctionMap(scale_by_time = true), dt = 1 / 4)
    sol2 = solve(prob_ode_linear, Euler(), dt = 1 / 4)
    @test sol[end] == sol2[end]
    @test sol(0.53) != sol2(0.53)
end

@testset "Array Discrete Problem" begin
    prob2 = DiscreteProblem(rand(4, 2), (0.0, 1.0))
    sol = solve(prob2, FunctionMap())
    @test sol[1] == sol[end]
    @test sol(0.5) == sol[1]

    @test_nowarn sol = solve(prob_ode_2Dlinear, FunctionMap())
    @test_nowarn sol2 = solve(prob_ode_2Dlinear, Euler(), dt = 1)
    sol = solve(prob_ode_2Dlinear, FunctionMap(scale_by_time = true), dt = 1 / 4)
    sol2 = solve(prob_ode_2Dlinear, Euler(), dt = 1 / 4)
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

@testset "DiscreteProblem in-place time definition" begin
    function f(un, p, t_np1)
        un + t_np1
    end

    prob2 = DiscreteProblem(f, 0, (0, 5))
    sol2 = solve(prob2, FunctionMap())

    @test sol2.u == [0, 1, 3, 6, 10, 15]

    function f!(u_np1, un, p, t_np1)
        u_np1[1] = un[1] + 1
        u_np1[2] = un[2] + t_np1
    end

    prob3 = DiscreteProblem(f!, [0, 0], (0, 5))
    sol3 = solve(prob3, FunctionMap())

    @test sol3.u == [[0, 0], [1, 1], [2, 3], [3, 6], [4, 10], [5, 15]]
end
