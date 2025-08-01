using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqCore
using Test
using Printf

"""
Allocation tests for OrdinaryDiffEqLowOrderRK solvers.
These tests verify that the step! operation does not allocate during stepping.
Based on testing: RK4 is allocation-free, but Euler needs fixed timestep.
"""

@testset "LowOrderRK Allocation Tests" begin
    # Test problem for adaptive methods
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    # Based on our testing, RK4 is allocation-free
    allocation_free_solvers = [RK4()]
    
    @testset "Known Allocation-Free LowOrderRK Solvers" begin
        for solver in allocation_free_solvers
            @testset "$(typeof(solver)) allocation test" begin
                integrator = init(prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                step!(integrator)  # Setup step may allocate
                
                # Test subsequent steps for zero allocations
                for i in 2:10
                    if integrator.t >= integrator.sol.prob.tspan[2]
                        break
                    end
                    alloc = @allocated step!(integrator)
                    @test alloc == 0
                end
            end
        end
    end
    
    # Test fixed timestep methods (like Euler) which require dt
    @testset "Fixed Timestep Methods" begin
        fixed_timestep_solvers = []
        
        try
            push!(fixed_timestep_solvers, Euler())
        catch
            # Euler may not be available
        end
        
        try
            push!(fixed_timestep_solvers, Midpoint())
        catch
            # Midpoint may not be available
        end
        
        try
            push!(fixed_timestep_solvers, Heun())
        catch
            # Heun may not be available
        end
        
        for solver in fixed_timestep_solvers
            @testset "$(typeof(solver)) fixed timestep allocation test" begin
                # Fixed timestep methods need dt specified
                integrator = init(prob, solver, dt=0.01, save_everystep=false, adaptive=false)
                step!(integrator)  # Setup step may allocate
                
                # Test subsequent steps for zero allocations
                for i in 2:10
                    if integrator.t >= integrator.sol.prob.tspan[2]
                        break
                    end
                    alloc = @allocated step!(integrator)
                    @test alloc == 0
                end
            end
        end
    end
end