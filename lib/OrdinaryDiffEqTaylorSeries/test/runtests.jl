using OrdinaryDiffEqTaylorSeries, ODEProblemLibrary, DiffEqDevTools
using SciMLBase
using Test
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @testset "Taylor2 Convergence Tests" begin
        # Test convergence
        dts = 2.0 .^ (-8:-4)
        testTol = 0.2
        sim = test_convergence(dts, prob_ode_linear, ExplicitTaylor2())
        @test sim.𝒪est[:final] ≈ 2 atol = testTol
        sim = test_convergence(dts, prob_ode_2Dlinear, ExplicitTaylor2())
        @test sim.𝒪est[:final] ≈ 2 atol = testTol
    end

    @testset "TaylorN Convergence Tests" begin
        # Test convergence
        dts = 2.0 .^ (-8:-4)
        testTol = 0.2
        for N in 3:4
            alg = ExplicitTaylor(order = Val(N))
            sim = test_convergence(dts, prob_ode_linear, alg)
            @test sim.𝒪est[:final] ≈ N atol = testTol
            sim = test_convergence(dts, prob_ode_2Dlinear, alg)
            @test sim.𝒪est[:final] ≈ N atol = testTol
        end
    end

    @testset "TaylorN Adaptive Tests" begin
        sol = solve(prob_ode_linear, ExplicitTaylor(order = Val(2)))
        @test length(sol) < 20
        @test SciMLBase.successful_retcode(sol)
    end

    # Test AutoSpecialize (default ODEProblem wraps in FunctionWrappers)
    # and FullSpecialize paths for IIP problems
    @testset "AutoSpecialize / FullSpecialize IIP" begin
        # IIP array problem
        function f_iip!(du, u, p, t)
            du[1] = -u[2]
            du[2] = u[1]
            return nothing
        end
        u0 = [1.0, 0.0]
        tspan = (0.0, 1.0)

        # AutoSpecialize (default) - uses FunctionWrappers, unwrapped_f needed
        prob_auto = ODEProblem(f_iip!, u0, tspan)
        # FullSpecialize - no wrapping
        prob_full = ODEProblem{true, SciMLBase.FullSpecialize}(f_iip!, u0, tspan)

        for prob in (prob_auto, prob_full)
            sol2 = solve(prob, ExplicitTaylor2(), dt = 0.01)
            @test SciMLBase.successful_retcode(sol2)
            @test length(sol2) > 1

            sol8 = solve(
                prob, ExplicitTaylor(order = Val(8)),
                abstol = 1.0e-12, reltol = 1.0e-12
            )
            @test SciMLBase.successful_retcode(sol8)
            @test length(sol8) > 1
        end

        # Verify both give similar results
        sol_auto = solve(
            prob_auto, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-12, reltol = 1.0e-12
        )
        sol_full = solve(
            prob_full, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-12, reltol = 1.0e-12
        )
        @test sol_auto.u[end] ≈ sol_full.u[end] atol = 1.0e-10
    end

    # Test OOP (out-of-place) with array state
    @testset "OOP Array Problems" begin
        function f_oop(u, p, t)
            return [-u[2], u[1]]
        end
        u0 = [1.0, 0.0]
        tspan = (0.0, 1.0)

        prob_oop = ODEProblem(f_oop, u0, tspan)

        sol2 = solve(prob_oop, ExplicitTaylor2(), dt = 0.01)
        @test SciMLBase.successful_retcode(sol2)
        @test length(sol2) > 1

        sol8 = solve(
            prob_oop, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-12, reltol = 1.0e-12
        )
        @test SciMLBase.successful_retcode(sol8)
        @test length(sol8) > 1

        # Check solution accuracy (harmonic oscillator: u1=cos(t), u2=sin(t))
        @test sol8.u[end][1] ≈ cos(1.0) atol = 1.0e-10
        @test sol8.u[end][2] ≈ sin(1.0) atol = 1.0e-10
    end

    # Dense output interpolation tests
    @testset "Taylor2 Dense Output" begin
        # Scalar OOP: u' = -u, u(0)=1 => u(t) = exp(-t)
        prob_scalar = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
        sol = solve(prob_scalar, ExplicitTaylor2(), dt = 0.01, dense = true)
        @test SciMLBase.successful_retcode(sol)
        # Check interpolation at intermediate points
        for t in 0.1:0.1:0.9
            @test sol(t) ≈ exp(-t) atol = 1.0e-3
        end

        # Array OOP: harmonic oscillator u1'=-u2, u2'=u1
        prob_arr = ODEProblem((u, p, t) -> [-u[2], u[1]], [1.0, 0.0], (0.0, 1.0))
        sol_arr = solve(prob_arr, ExplicitTaylor2(), dt = 0.01, dense = true)
        @test SciMLBase.successful_retcode(sol_arr)
        t_mid = 0.5
        @test sol_arr(t_mid)[1] ≈ cos(t_mid) atol = 1.0e-3
        @test sol_arr(t_mid)[2] ≈ sin(t_mid) atol = 1.0e-3

        # idxs parameter
        @test sol_arr(t_mid, idxs = 1) ≈ cos(t_mid) atol = 1.0e-3
        @test sol_arr(t_mid, idxs = 2) ≈ sin(t_mid) atol = 1.0e-3

        # IIP: same harmonic oscillator
        function f_dense_iip!(du, u, p, t)
            du[1] = -u[2]
            du[2] = u[1]
            return nothing
        end
        prob_iip = ODEProblem{true, SciMLBase.FullSpecialize}(
            f_dense_iip!, [1.0, 0.0], (0.0, 1.0)
        )
        sol_iip = solve(prob_iip, ExplicitTaylor2(), dt = 0.01, dense = true)
        @test SciMLBase.successful_retcode(sol_iip)
        @test sol_iip(t_mid)[1] ≈ cos(t_mid) atol = 1.0e-3
        @test sol_iip(t_mid)[2] ≈ sin(t_mid) atol = 1.0e-3
    end

    @testset "Taylor2 Dense Convergence" begin
        dts = 2.0 .^ (-8:-4)
        testTol = 0.2
        sim = test_convergence(
            dts, prob_ode_linear, ExplicitTaylor2(), dense_errors = true
        )
        @test sim.𝒪est[:L2] ≈ 2 atol = testTol
        sim = test_convergence(
            dts, prob_ode_2Dlinear, ExplicitTaylor2(), dense_errors = true
        )
        @test sim.𝒪est[:L2] ≈ 2 atol = testTol
    end

    @testset "TaylorN Dense Output" begin
        # Scalar OOP: u' = -u, u(0)=1 => u(t) = exp(-t)
        prob_scalar = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))

        for N in (4, 8)
            alg = ExplicitTaylor(order = Val(N))
            sol = solve(prob_scalar, alg, abstol = 1.0e-12, reltol = 1.0e-12, dense = true)
            @test SciMLBase.successful_retcode(sol)
            # High-order Taylor should give very accurate interpolation
            tol = N >= 8 ? 1.0e-10 : 1.0e-6
            for t in 0.1:0.1:0.9
                @test sol(t) ≈ exp(-t) atol = tol
            end
        end

        # Array OOP: harmonic oscillator
        prob_arr = ODEProblem((u, p, t) -> [-u[2], u[1]], [1.0, 0.0], (0.0, 1.0))
        sol8 = solve(
            prob_arr, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-12, reltol = 1.0e-12, dense = true
        )
        @test SciMLBase.successful_retcode(sol8)
        t_mid = 0.5
        @test sol8(t_mid)[1] ≈ cos(t_mid) atol = 1.0e-10
        @test sol8(t_mid)[2] ≈ sin(t_mid) atol = 1.0e-10

        # idxs parameter
        @test sol8(t_mid, idxs = 1) ≈ cos(t_mid) atol = 1.0e-10
        @test sol8(t_mid, idxs = 2) ≈ sin(t_mid) atol = 1.0e-10

        # IIP: harmonic oscillator with AutoSpecialize
        function f_taylor_dense_iip!(du, u, p, t)
            du[1] = -u[2]
            du[2] = u[1]
            return nothing
        end
        prob_iip = ODEProblem(f_taylor_dense_iip!, [1.0, 0.0], (0.0, 1.0))
        sol_iip = solve(
            prob_iip, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-12, reltol = 1.0e-12, dense = true
        )
        @test SciMLBase.successful_retcode(sol_iip)
        @test sol_iip(t_mid)[1] ≈ cos(t_mid) atol = 1.0e-10
        @test sol_iip(t_mid)[2] ≈ sin(t_mid) atol = 1.0e-10

        # idxs on IIP
        @test sol_iip(t_mid, idxs = 1) ≈ cos(t_mid) atol = 1.0e-10
        @test sol_iip(t_mid, idxs = 2) ≈ sin(t_mid) atol = 1.0e-10
    end

    @testset "TaylorN Dense Convergence" begin
        dts = 2.0 .^ (-8:-4)
        testTol = 0.3
        for N in 3:4
            alg = ExplicitTaylor(order = Val(N))
            sim = test_convergence(
                dts, prob_ode_linear, alg, dense_errors = true
            )
            @test sim.𝒪est[:L2] ≈ N atol = testTol
            sim = test_convergence(
                dts, prob_ode_2Dlinear, alg, dense_errors = true
            )
            @test sim.𝒪est[:L2] ≈ N atol = testTol
        end
    end

    # Test IIP with a nonlinear system (Henon-Heiles-style)
    @testset "Nonlinear IIP System" begin
        function henon_heiles!(du, u, p, t)
            du[1] = -u[3] * (1 + 2u[4])
            du[2] = -u[4] - (u[3]^2 - u[4]^2)
            du[3] = u[1]
            du[4] = u[2]
            return nothing
        end
        u0 = [0.0, 0.5, 0.1, 0.0]
        tspan = (0.0, 10.0)

        prob = ODEProblem{true, SciMLBase.FullSpecialize}(henon_heiles!, u0, tspan)
        sol = solve(
            prob, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-14, reltol = 1.0e-14
        )
        @test SciMLBase.successful_retcode(sol)
        @test length(sol) > 1

        # Also test with AutoSpecialize
        prob_auto = ODEProblem(henon_heiles!, u0, tspan)
        sol_auto = solve(
            prob_auto, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-14, reltol = 1.0e-14
        )
        @test SciMLBase.successful_retcode(sol_auto)
        @test sol_auto.u[end] ≈ sol.u[end] atol = 1.0e-10
    end
end

# Run QA tests (JET, Aqua)
if TEST_GROUP != "FUNCTIONAL" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
