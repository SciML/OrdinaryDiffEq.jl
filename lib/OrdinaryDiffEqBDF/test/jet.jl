using Pkg
Pkg.add("JET")

import OrdinaryDiffEqBDF
using OrdinaryDiffEqBDF
using OrdinaryDiffEqCore
using DiffEqBase: SplitODEProblem, DAEProblem
using JET
using Test

@testset "JET Tests" begin
    # Test package for typos - now passing
    test_package(
        OrdinaryDiffEqBDF, target_modules = (OrdinaryDiffEqBDF,), mode = :typo
    )

    # Test individual solver type stability
    @testset "Solver Type Stability Tests" begin
        # Test problem - use a simple linear problem for stiff solvers
        linear_prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))

        # Split problem for SBDF solvers (which require SplitODEProblem)
        split_prob = SplitODEProblem((u, p, t) -> -u, (u, p, t) -> 0.0, 1.0, (0.0, 1.0))

        # DAE problem for DAE solvers
        function simple_dae!(resid, du, u, p, t)
            resid[1] = -u[1] - du[1]
        end
        u0 = [1.0]
        du0 = [-1.0]
        dae_prob = DAEProblem(simple_dae!, du0, u0, (0.0, 1.0))

        # Regular BDF solvers (ODEProblem)
        regular_bdf_solvers = [
            ABDF2(), QNDF1(), QBDF1(), QNDF2(), QBDF2(), QNDF(), QBDF(), FBDF(),
            MEBDF2(),
        ]
        # After the controller-cache refactor every regular BDF init is
        # type-stable for `init`.
        stable_bdf_solvers = regular_bdf_solvers

        # DAE solvers (DAEProblem)
        dae_solvers = [DABDF2(), DImplicitEuler(), DFBDF()]

        # Test SBDF solvers separately with required order parameter and SplitODEProblem
        sbdf_solvers = [IMEXEuler(), IMEXEulerARK(), SBDF(order = 2), SBDF(order = 3), SBDF(order = 4), SBDF2(), SBDF3(), SBDF4()]

        for solver in regular_bdf_solvers
            @testset "$(typeof(solver)) type stability" begin
                sbroken = solver ∉ stable_bdf_solvers
                @test_opt broken = sbroken init(linear_prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                integrator = init(linear_prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                @test_opt broken = true step!(integrator)
            end
        end

        for solver in dae_solvers
            @testset "$(typeof(solver)) DAE type stability" begin
                @test_opt broken = true init(dae_prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                integrator = init(dae_prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                @test_opt broken = true step!(integrator)
            end
        end

        for solver in sbdf_solvers
            @testset "$(typeof(solver)) type stability" begin
                @test_opt init(split_prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                integrator = init(split_prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                @test_opt broken = true step!(integrator)
            end
        end
    end
end
