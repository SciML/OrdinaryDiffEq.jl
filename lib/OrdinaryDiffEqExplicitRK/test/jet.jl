import OrdinaryDiffEqExplicitRK
using OrdinaryDiffEqExplicitRK
using OrdinaryDiffEqCore
using JET
using Test

@testset "JET Tests" begin
    # Test package for typos
    test_package(
        OrdinaryDiffEqExplicitRK, target_defined_modules = true, mode = :typo)
    
    # Test individual solver type stability
    @testset "Solver Type Stability Tests" begin
        # Test problem
        function simple_system!(du, u, p, t)
            du[1] = -0.5 * u[1]
            du[2] = -1.5 * u[2]
        end
        prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
        
        # Test all exported ExplicitRK solvers
        explicit_rk_solvers = [ExplicitRK()]
        
        for solver in explicit_rk_solvers
            @testset "$(typeof(solver)) type stability" begin
                @test_opt init(prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                integrator = init(prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
                @test_opt step!(integrator)
            end
        end
    end
end
