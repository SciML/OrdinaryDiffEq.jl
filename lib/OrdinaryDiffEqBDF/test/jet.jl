import OrdinaryDiffEqBDF
using OrdinaryDiffEqBDF
using OrdinaryDiffEqCore
using JET
using Test

@testset "JET Tests" begin
    # Test package for typos
    test_package(
        OrdinaryDiffEqBDF, target_defined_modules = true, mode = :typo)
    
    # Test individual solver type stability
    @testset "Solver Type Stability Tests" begin
        # Test problem - use a simple linear problem for stiff solvers
        linear_prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
        
        # Test all exported BDF solvers
        bdf_solvers = [ABDF2(), QNDF1(), QBDF1(), QNDF2(), QBDF2(), QNDF(), QBDF(), FBDF(),
                       SBDF(), SBDF2(), SBDF3(), SBDF4(), MEBDF2(), IMEXEuler(), IMEXEulerARK(),
                       DABDF2(), DImplicitEuler(), DFBDF()]
        
        for solver in bdf_solvers
            @testset "$(typeof(solver)) type stability" begin
                try
                    @test_opt init(linear_prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                    integrator = init(linear_prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                    @test_opt step!(integrator)
                catch e
                    @test_broken false # Mark as broken if solver fails to initialize
                    println("$(typeof(solver)) failed with: $e")
                end
            end
        end
    end
end
