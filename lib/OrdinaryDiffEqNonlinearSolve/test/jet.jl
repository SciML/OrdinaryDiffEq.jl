using Pkg
Pkg.add("JET")

import OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqBDF
using LinearAlgebra: Diagonal
using ADTypes
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqNonlinearSolve, target_modules = (OrdinaryDiffEqNonlinearSolve,), mode = :typo
    )

    # Test that DAE initialization is type-stable for a sampling of algorithms
    @testset "Solver Type Stability Tests" begin

        # Mass matrix problem
        function dae(du, u, p, t)
            du[1] = u[2]
            du[2] = exp(u[1]) + exp(u[2]) # deliberately not satisfiable
            return nothing
        end
        function dae_oop(u, p, t)
            return [u[2], exp(u[1]) + exp(u[2])]
        end
        dae_func = ODEFunction(dae, mass_matrix = Diagonal([1.0, 0.0]))
        dae_oop_func = ODEFunction(dae_oop, mass_matrix = Diagonal([1.0, 0.0]))
        mm_prob = ODEProblem(dae_func, [0.0, 1.0], (0.0, 1.0))
        mm_oop_prob = ODEProblem(dae_oop_func, [0.0, 1.0], (0.0, 1.0))
        ad = AutoForwardDiff(chunksize = 2)

        # DAE problem for DAE solvers
        function simple_dae!(resid, du, u, p, t)
            resid[1] = -u[1] - du[1]
        end
        u0 = [1.0]
        du0 = [-1.0]
        dae_prob = DAEProblem(simple_dae!, du0, u0, (0.0, 1.0))


        # Solvers which are not fully type-stable for initialization
        iip_unstable_ode_solvers = [
            ABDF2(), QNDF1(), QBDF1(), QNDF2(), QBDF2(), QNDF(), QBDF(), FBDF(), MEBDF2(),
            Rosenbrock23(), Rosenbrock32(), Rodas4(), Rodas4P(), Rodas5(), Rodas5P(),
        ]
        # Some of these are type-stable for the initialization step, but not all
        iip_stable_ode_solvers = [
            Rosenbrock23(autodiff = ad), Rosenbrock32(autodiff = ad), Rodas4(autodiff = ad),
            Rodas4P(autodiff = ad), Rodas5(autodiff = ad), Rodas5P(autodiff = ad),
        ]

        # Mass matrix ODE problem
        @testset "iip_ode_solvers" begin
            for solver in iip_unstable_ode_solvers
                @testset "$(typeof(solver)) type stability" begin
                    @test_opt broken = true init(mm_prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                end
            end
            for solver in iip_stable_ode_solvers
                @testset "$(typeof(solver)) type stability" begin
                    @test_opt init(mm_prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                end
            end
        end

        # Solvers which are not fully type-stable for initialization
        oop_unstable_ode_solvers = [
            ABDF2(), QNDF1(), QBDF1(), QNDF2(), QBDF2(), QNDF(), QBDF(), FBDF(), MEBDF2(),
            Rosenbrock23(), Rosenbrock32(), Rodas4(), Rodas4P(), Rodas5(), Rodas5P(),
        ]
        # Some of these are type-stable for the initialization step, but not all
        oop_stable_ode_solvers = [
            Rosenbrock23(autodiff = ad), Rosenbrock32(autodiff = ad), Rodas4(autodiff = ad),
            Rodas4P(autodiff = ad), Rodas5(autodiff = ad), Rodas5P(autodiff = ad),
        ]

        # Mass matrix ODE problem
        @testset "oop_ode_solvers" begin
            for solver in oop_unstable_ode_solvers
                @testset "$(typeof(solver)) type stability" begin
                    @test_opt broken = true init(mm_prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                end
            end
            for solver in oop_stable_ode_solvers
                @testset "$(typeof(solver)) type stability" begin
                    @test_opt init(mm_prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                end
            end
        end

        # DAE solvers (DAEProblem)
        dae_solvers = [DABDF2(), DImplicitEuler(), DFBDF()]


        for solver in dae_solvers
            @testset "$(typeof(solver)) DAE type stability" begin
                @test_opt broken = true init(dae_prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
            end
        end

    end
end
