using OrdinaryDiffEqVerner
using OrdinaryDiffEqCore
using Test
using Printf

"""
Allocation tests for OrdinaryDiffEqVerner solvers.
These tests verify that the step! operation does not allocate during stepping.
"""

@testset "Verner Allocation Tests" begin
    # Test problem
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    # Based on our testing, Verner solvers are allocation-free
    # Note: Vern6-Vern9 are in HighOrderRK, so test others here
    available_solvers = []
    
    # Try to add available Verner solvers
    try
        push!(available_solvers, Vern6())
    catch
        # Vern6 may be in different package
    end
    
    try
        push!(available_solvers, Vern7())
    catch
        # Vern7 may be in different package
    end
    
    if !isempty(available_solvers)
        @testset "Known Allocation-Free Verner Solvers" begin
            for solver in available_solvers
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
    else
        @test_skip "No Verner solvers available in this package"
    end
end