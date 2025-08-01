using OrdinaryDiffEqSSPRK
using OrdinaryDiffEqCore
using Test
using Printf

"""
Allocation tests for OrdinaryDiffEqSSPRK solvers.
These tests verify that the step! operation does not allocate during stepping.
Based on testing: SSPRK43 is allocation-free, but SSPRK22 and SSPRK33 are allocating.
"""

@testset "SSPRK Allocation Tests" begin
    # Test problem
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    # Based on our testing, SSPRK43 is allocation-free
    allocation_free_solvers = []
    try
        push!(allocation_free_solvers, SSPRK43())
    catch
        # SSPRK43 may not be available
    end
    
    @testset "Known Allocation-Free SSPRK Solvers" begin
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
    
    # Test currently allocating SSPRK solvers with @test_broken
    allocating_solvers = []
    try
        push!(allocating_solvers, SSPRK22())
    catch
        # SSPRK22 may not be available
    end
    
    try  
        push!(allocating_solvers, SSPRK33())
    catch
        # SSPRK33 may not be available
    end
    
    @testset "Currently Allocating SSPRK Solvers (@test_broken)" begin
        for solver in allocating_solvers
            @testset "$(typeof(solver)) allocation test (broken)" begin
                integrator = init(prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                step!(integrator)  # Setup step may allocate
                
                # These tests are expected to fail until allocation issues are resolved
                for i in 2:6
                    if integrator.t >= integrator.sol.prob.tspan[2]
                        break
                    end
                    alloc = @allocated step!(integrator)
                    @test_broken alloc == 0
                end
            end
        end
    end
end