import OrdinaryDiffEqFunctionMap
import OrdinaryDiffEqCore
using Test
using SciMLBase
using SciMLBase: solve, init, DiscreteProblem

const FunctionMap = OrdinaryDiffEqFunctionMap.FunctionMap

# Helper functions to check algorithm properties regardless of module context
is_functionmap(alg) = typeof(alg).name.name == :FunctionMap
function get_scale_by_time(alg)
    # Access the type parameter directly since it's FunctionMap{scale_by_time}
    T = typeof(alg)
    # The parameter is stored as a type parameter
    return T.parameters[1]
end

@testset "DiscreteProblem Default Algorithm" begin
    # Test scalar DiscreteProblem
    f(u, p, t) = 1.01 * u
    prob_scalar = DiscreteProblem(f, 0.5, (0.0, 1.0))

    @testset "Scalar DiscreteProblem" begin
        # Test solve without explicit algorithm
        sol = solve(prob_scalar)
        @test is_functionmap(sol.alg)
        @test get_scale_by_time(sol.alg) == false
        @test length(sol.u) > 1

        # Test init without explicit algorithm
        integrator = init(prob_scalar)
        @test is_functionmap(integrator.alg)
        @test get_scale_by_time(integrator.alg) == false
    end

    # Test array DiscreteProblem
    function f_array!(du, u, p, t)
        du[1] = 1.01 * u[1]
        du[2] = 0.99 * u[2]
    end
    prob_array = DiscreteProblem(f_array!, [0.5, 1.0], (0.0, 1.0))

    @testset "Array DiscreteProblem" begin
        # Test solve without explicit algorithm
        sol = solve(prob_array)
        @test is_functionmap(sol.alg)
        @test get_scale_by_time(sol.alg) == false
        @test length(sol.u) > 1

        # Test init without explicit algorithm
        integrator = init(prob_array)
        @test is_functionmap(integrator.alg)
        @test get_scale_by_time(integrator.alg) == false
    end

    # Test that explicit algorithm specification still works
    @testset "Explicit FunctionMap specification" begin
        sol1 = solve(prob_scalar, FunctionMap())
        @test is_functionmap(sol1.alg)
        @test get_scale_by_time(sol1.alg) == false

        sol2 = solve(prob_scalar, FunctionMap(scale_by_time = true), dt = 0.1)
        @test is_functionmap(sol2.alg)
        @test get_scale_by_time(sol2.alg) == true

        integrator1 = init(prob_scalar, FunctionMap())
        @test is_functionmap(integrator1.alg)
        @test get_scale_by_time(integrator1.alg) == false

        integrator2 = init(prob_scalar, FunctionMap(scale_by_time = true), dt = 0.1)
        @test is_functionmap(integrator2.alg)
        @test get_scale_by_time(integrator2.alg) == true
    end

    # Test that the default behaves correctly with different problem types
    @testset "DiscreteProblem with integer time" begin
        henon_map!(u_next, u, p, t) = begin
            u_next[1] = 1 + u[2] - 1.4 * u[1]^2
            u_next[2] = 0.3 * u[1]
        end

        prob_int = DiscreteProblem(henon_map!, [0.5, 0.5], (0, 10))

        sol = solve(prob_int)
        @test is_functionmap(sol.alg)
        @test eltype(sol.t) <: Integer

        integrator = init(prob_int)
        @test is_functionmap(integrator.alg)
    end
end
