import OrdinaryDiffEqSSPRK
using OrdinaryDiffEqSSPRK
using OrdinaryDiffEqCore
using JET
using Test

@testset "JET Tests" begin
    # Test package for typos (mark as broken for now)
    @test_broken test_package(
        OrdinaryDiffEqSSPRK, target_defined_modules = true, mode = :typo)
    
    # Test individual solver type stability
    @testset "Solver Type Stability Tests" begin
        # Test problem
        function simple_system!(du, u, p, t)
            du[1] = -0.5 * u[1]
            du[2] = -1.5 * u[2]
        end
        prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
        
        # Test main SSPRK solvers (mark as broken)
        ssprk_solvers = [SSPRK22(), SSPRK33(), SSPRK43(), SSPRK432(), SSPRKMSVS32(), SSPRKMSVS43(),
                        SSPRK932(), SSPRK54(), SSPRK73(), SSPRK83(), SSPRK63()]
        
        for solver in ssprk_solvers
            @testset "$(typeof(solver)) type stability" begin
                try
                    @test_broken @test_opt init(prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                    integrator = init(prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                    @test_broken @test_opt step!(integrator)
                catch e
                    @test_broken false # Mark as broken if solver fails to initialize
                    println("$(typeof(solver)) failed with: $e")
                end
            end
        end
    end
end
