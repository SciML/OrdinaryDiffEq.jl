using OrdinaryDiffEqFunctionMap
using OrdinaryDiffEqCore
using Test
using SciMLBase

@testset "DiscreteProblem Default Algorithm" begin
    # Test scalar DiscreteProblem
    f(u, p, t) = 1.01 * u
    prob_scalar = DiscreteProblem(f, 0.5, (0.0, 1.0))
    
    @testset "Scalar DiscreteProblem" begin
        # Test solve without explicit algorithm
        sol = solve(prob_scalar)
        @test typeof(sol.alg).name.name == :FunctionMap
        @test sol.alg == FunctionMap(scale_by_time=false)
        @test length(sol.u) > 1
        
        # Test init without explicit algorithm
        integrator = init(prob_scalar)
        @test typeof(integrator.alg).name.name == :FunctionMap
        @test integrator.alg == FunctionMap(scale_by_time=false)
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
        @test typeof(sol.alg).name.name == :FunctionMap
        @test sol.alg == FunctionMap(scale_by_time=false)
        @test length(sol.u) > 1
        
        # Test init without explicit algorithm
        integrator = init(prob_array)
        @test typeof(integrator.alg).name.name == :FunctionMap
        @test integrator.alg == FunctionMap(scale_by_time=false)
    end
    
    # Test that explicit algorithm specification still works
    @testset "Explicit FunctionMap specification" begin
        sol1 = solve(prob_scalar, FunctionMap())
        @test sol1.alg == FunctionMap(scale_by_time=false)
        
        sol2 = solve(prob_scalar, FunctionMap(scale_by_time=true), dt=0.1)
        @test sol2.alg == FunctionMap(scale_by_time=true)
        
        integrator1 = init(prob_scalar, FunctionMap())
        @test integrator1.alg == FunctionMap(scale_by_time=false)
        
        integrator2 = init(prob_scalar, FunctionMap(scale_by_time=true), dt=0.1)
        @test integrator2.alg == FunctionMap(scale_by_time=true)
    end
    
    # Test that the default behaves correctly with different problem types
    @testset "DiscreteProblem with integer time" begin
        henon_map!(u_next, u, p, t) = begin
            u_next[1] = 1 + u[2] - 1.4 * u[1]^2
            u_next[2] = 0.3 * u[1]
        end
        
        prob_int = DiscreteProblem(henon_map!, [0.5, 0.5], (0, 10))
        
        sol = solve(prob_int)
        @test typeof(sol.alg).name.name == :FunctionMap
        @test eltype(sol.t) <: Integer
        
        integrator = init(prob_int)
        @test typeof(integrator.alg).name.name == :FunctionMap
    end
end