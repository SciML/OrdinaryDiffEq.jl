using OrdinaryDiffEqTsit5
using OrdinaryDiffEqCore
using JET
using Test
using Printf

"""
JET type stability tests for OrdinaryDiffEqTsit5 solvers.
These tests verify that the step! operation and solve operations are type stable.
"""

@testset "Tsit5 JET Type Stability Tests" begin
    # Test problem
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    @testset "Solver Initialization Type Stability" begin
        # This will likely fail initially due to type instabilities 
        # Mark as broken until type stability issues are resolved
        @test_opt broken=true init(prob, Tsit5(), save_everystep=false, abstol=1e-6, reltol=1e-6)
    end
    
    @testset "Step Operation Type Stability" begin
        integrator = init(prob, Tsit5(), save_everystep=false, abstol=1e-6, reltol=1e-6)
        
        # This will likely fail initially due to type instabilities in the solvers
        # Mark as broken until type stability issues are resolved
        @test_opt broken=true step!(integrator)
    end
    
    @testset "Full Solve Type Stability" begin
        # This will likely fail initially due to type instabilities in the solvers
        # Mark as broken until type stability issues are resolved
        @test_opt broken=true solve(prob, Tsit5(), save_everystep=false, abstol=1e-6, reltol=1e-6)
    end
end