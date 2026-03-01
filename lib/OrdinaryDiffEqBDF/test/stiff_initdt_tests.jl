using Test
using OrdinaryDiffEqBDF
using OrdinaryDiffEqSDIRK

@testset "StiffInitDt Algorithm" begin

    @testset "Simple exponential decay (in-place)" begin
        function f_exp!(du, u, p, t)
            du[1] = -u[1]
        end
        prob = ODEProblem(f_exp!, [1.0], (0.0, 10.0))
        sol = solve(prob, FBDF())
        @test sol.retcode == ReturnCode.Success
        @test isapprox(sol.u[end][1], exp(-10.0), rtol = 1.0e-2)
    end

    @testset "Simple exponential decay (out-of-place)" begin
        f_exp(u, p, t) = [-u[1]]
        prob = ODEProblem(f_exp, [1.0], (0.0, 10.0))
        sol = solve(prob, FBDF())
        @test sol.retcode == ReturnCode.Success
        @test isapprox(sol.u[end][1], exp(-10.0), rtol = 1.0e-2)
    end

    @testset "Multi-scale with zero states and tiny abstol (in-place)" begin
        # This is the pathological case from issue #1496:
        # Some species start at 0 with tiny abstol, causing the Hairer algorithm
        # to produce catastrophically small initial dt.
        function f_ms!(du, u, p, t)
            du[1] = -100.0 * u[1] + 1.0
            du[2] = 0.1 * u[1]
            du[3] = 1.84  # mimics TH2S
        end
        u0 = [1.0, 0.0, 0.0]
        prob = ODEProblem(f_ms!, u0, (0.0, 100.0))

        # With FBDF (StiffInitDt) - should succeed efficiently
        sol_fbdf = solve(prob, FBDF(), abstol = 1.0e-15, reltol = 1.0e-8)
        @test sol_fbdf.retcode == ReturnCode.Success
        @test sol_fbdf.t[end] == 100.0

        # The StiffInitDt should produce an initial dt much larger than the
        # catastrophically small ~5e-25 that the Hairer algorithm would give.
        # With abstol=1e-15 and u0=0, the StiffInitDt computes ~3.5e-15
        # (10 orders of magnitude larger).
        dt_fbdf = sol_fbdf.t[2] - sol_fbdf.t[1]
        @test dt_fbdf > 1.0e-20  # Much larger than the ~5e-25 Hairer would give
    end

    @testset "Multi-scale with zero states and tiny abstol (out-of-place)" begin
        function f_ms(u, p, t)
            return [-100.0 * u[1] + 1.0, 0.1 * u[1], 1.84]
        end
        u0 = [1.0, 0.0, 0.0]
        prob = ODEProblem(f_ms, u0, (0.0, 100.0))

        sol = solve(prob, FBDF(), abstol = 1.0e-15, reltol = 1.0e-8)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] == 100.0
    end

    @testset "Stiff Robertson problem" begin
        # Classic stiff test problem
        function robertson!(du, u, p, t)
            du[1] = -0.04 * u[1] + 1.0e4 * u[2] * u[3]
            du[2] = 0.04 * u[1] - 1.0e4 * u[2] * u[3] - 3.0e7 * u[2]^2
            du[3] = 3.0e7 * u[2]^2
        end
        u0 = [1.0, 0.0, 0.0]
        prob = ODEProblem(robertson!, u0, (0.0, 1.0e5))

        sol = solve(prob, FBDF(), abstol = 1.0e-8, reltol = 1.0e-8)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] == 1.0e5

        # Conservation law: u1 + u2 + u3 = 1
        @test isapprox(sum(sol.u[end]), 1.0, atol = 1.0e-6)
    end

    @testset "Backward integration" begin
        function f_back!(du, u, p, t)
            du[1] = -u[1]
        end
        prob = ODEProblem(f_back!, [exp(-10.0)], (10.0, 0.0))
        sol = solve(prob, FBDF())
        @test sol.retcode == ReturnCode.Success
        @test isapprox(sol.u[end][1], 1.0, rtol = 2.0e-2)
    end

    @testset "User-specified dt still works" begin
        function f_dt!(du, u, p, t)
            du[1] = -u[1]
        end
        prob = ODEProblem(f_dt!, [1.0], (0.0, 1.0))
        # When user specifies dt, StiffInitDt should not be called
        sol = solve(prob, FBDF(), dt = 0.01)
        @test sol.retcode == ReturnCode.Success
    end

    @testset "QNDF with multi-scale problem" begin
        function f_qndf!(du, u, p, t)
            du[1] = -100.0 * u[1] + 1.0
            du[2] = 0.0  # starts at 0, stays at 0
        end
        u0 = [1.0, 0.0]
        prob = ODEProblem(f_qndf!, u0, (0.0, 1.0))

        sol = solve(prob, QNDF(), abstol = 1.0e-12, reltol = 1.0e-8)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] == 1.0
    end

    @testset "ImplicitEuler with StiffInitDt" begin
        # ImplicitEuler now also uses StiffInitDt since all implicit methods do
        function f_ie!(du, u, p, t)
            du[1] = -100.0 * u[1] + 1.0
            du[2] = 0.0
        end
        u0 = [1.0, 0.0]
        prob = ODEProblem(f_ie!, u0, (0.0, 1.0))

        sol = solve(prob, ImplicitEuler(), abstol = 1.0e-12, reltol = 1.0e-8)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] == 1.0
    end
end
