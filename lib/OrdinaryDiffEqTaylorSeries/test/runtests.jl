using SciMLTesting
using OrdinaryDiffEqTaylorSeries, ODEProblemLibrary, DiffEqDevTools
using SciMLBase
using Test
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @testset "ExplicitTaylor2 Convergence Tests" begin
        # Test convergence
        dts = 2.0 .^ (-8:-4)
        testTol = 0.2
        sim = test_convergence(dts, prob_ode_linear, ExplicitTaylor2())
        @test sim.𝒪est[:final] ≈ 2 atol = testTol
        sim = test_convergence(dts, prob_ode_2Dlinear, ExplicitTaylor2())
        @test sim.𝒪est[:final] ≈ 2 atol = testTol
    end

    @testset "ExplicitTaylorN Convergence Tests" begin
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

    @testset "ExplicitTaylorAdaptiveOrder Tests" begin
        sol = solve(prob_ode_linear, ExplicitTaylorAdaptiveOrder(min_order = Val(6), max_order = Val(10)))
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)
        # a step count and a retcode alone are satisfied by a solver that never
        # leaves `u0`, which is exactly what issue #3976 was
        @test sol.u[end] ≈ prob_ode_linear.f.analytic(
            prob_ode_linear.u0, prob_ode_linear.p, sol.t[end]
        ) atol = 1.0e-6
    end

    # Issue #3976. Every order's jet writes into one Taylor buffer sized at
    # `max_order`. Storing a lower-order `TaylorScalar` there went through
    # TaylorDiff's number-promotion constructor, which keeps the constant term and
    # drops every derivative, so the step collapsed to `u = uprev`, the error
    # estimate was exactly zero and the solver reported Success after a handful of
    # enormous steps.
    @testset "ExplicitTaylorAdaptiveOrder keeps every Taylor coefficient" begin
        TaylorDiff = OrdinaryDiffEqTaylorSeries.TaylorDiff
        function f_exp!(du, u, p, t)
            du .= u
            return nothing
        end
        prob = ODEProblem{true, SciMLBase.FullSpecialize}(f_exp!, [1.0, 2.0], (0.0, 1.0))
        integ = init(
            prob, ExplicitTaylorAdaptiveOrder(min_order = Val(2), max_order = Val(6)),
            abstol = 1.0e-10, reltol = 1.0e-10
        )
        cache = integ.cache
        # `jets[i]` has order `min_order + i - 1`, and for u' = u it owes the buffer
        # that many nonzero coefficients on top of the value
        for (i, jet) in enumerate(cache.jets)
            fill!(cache.utaylor, zero(eltype(cache.utaylor)))
            jet(cache.utaylor, cache.coeffs[i], integ.uprev, integ.t)
            @test count(!iszero, TaylorDiff.flatten(cache.utaylor[1])) == i + 2
        end
    end

    @testset "ExplicitTaylorAdaptiveOrder accuracy tracks tolerance" begin
        for prob in (prob_ode_linear, prob_ode_2Dlinear)
            for tol in (1.0e-8, 1.0e-10)
                sol = solve(prob, ExplicitTaylorAdaptiveOrder(), abstol = tol, reltol = tol)
                @test SciMLBase.successful_retcode(sol)
                @test sol.u[end] != prob.u0
                exact = prob.f.analytic(prob.u0, prob.p, sol.t[end])
                @test maximum(abs.(sol.u[end] .- exact)) < 1.0e-6 * maximum(abs.(exact))
            end
        end

        # out-of-place array problems additionally hit an UndefRefError while
        # saving: `initialize!` sized `integrator.k` but `perform_step!` never
        # filled it
        prob_arr = ODEProblem((u, p, t) -> [-u[2], u[1]], [1.0, 0.0], (0.0, 1.0))
        sol = solve(
            prob_arr, ExplicitTaylorAdaptiveOrder(),
            abstol = 1.0e-10, reltol = 1.0e-10
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end] ≈ [cos(1.0), sin(1.0)] atol = 1.0e-8
        @test sol(0.5) ≈ [cos(0.5), sin(0.5)] atol = 1.0e-8

        # the order window has to leave room for the extra jet the error estimate
        # rides on, otherwise the order loop never runs and a stale estimate is
        # what the step is accepted on
        @test_throws ArgumentError solve(
            prob_ode_linear,
            ExplicitTaylorAdaptiveOrder(min_order = Val(4), max_order = Val(4))
        )
    end

    # `test_convergence` solves with `adaptive = false`, which pins the order at
    # `max_order - 1` for the whole run. The jet order then equals the buffer
    # order, so this covers the top of the window only.
    @testset "ExplicitTaylorAdaptiveOrder Convergence Tests" begin
        dts = 2.0 .^ (-6:-2)
        testTol = 0.2
        for N in 3:5
            alg = ExplicitTaylorAdaptiveOrder(min_order = Val(N - 2), max_order = Val(N))
            sim = test_convergence(dts, prob_ode_linear, alg)
            @test sim.𝒪est[:final] ≈ N atol = testTol
            sim = test_convergence(dts, prob_ode_2Dlinear, alg)
            @test sim.𝒪est[:final] ≈ N atol = testTol
        end
    end

    # The mismatch of #3976 only exists below the top of the window, so hold the
    # order there by hand and check that the method still converges at the order
    # its jet claims. With the coefficients dropped the step is a no-op and the
    # error stops depending on `dt` at all.
    @testset "ExplicitTaylorAdaptiveOrder Convergence Below max_order" begin
        alg = ExplicitTaylorAdaptiveOrder(min_order = Val(2), max_order = Val(6))
        dts = 2.0 .^ (-8:-4)
        for prob in (prob_ode_linear, prob_ode_2Dlinear)
            exact = prob.f.analytic(prob.u0, prob.p, prob.tspan[end])
            errors = map(dts) do dt
                integrator = init(prob, alg; dt, adaptive = false)
                integrator.cache.current_order[] = 3
                solve!(integrator)
                maximum(abs.(integrator.sol.u[end] .- exact))
            end
            slopes = diff(log2.(errors)) ./ diff(log2.(dts))
            @test sum(slopes) / length(slopes) ≈ 4 atol = 0.2
        end
    end

    # Dense output on the in-place cache reads `max_order` coefficients back out
    # of `integrator.k`, which the adaptive-order step has to fill and the
    # interpolant has to read at the buffer order rather than at `min_order`.
    @testset "ExplicitTaylorAdaptiveOrder Dense Output" begin
        for prob in (prob_ode_linear, prob_ode_2Dlinear)
            sol = solve(
                prob, ExplicitTaylorAdaptiveOrder(),
                abstol = 1.0e-10, reltol = 1.0e-10, dense = true
            )
            for t in (0.25, 0.5, 0.75)
                exact = prob.f.analytic(prob.u0, prob.p, t)
                @test maximum(abs.(sol(t) .- exact)) < 1.0e-8 * maximum(abs.(exact))
            end
        end
    end

    # Issue #3976: once the order dropped off `max_order` the step froze, so the
    # solver walked to the end of the span in a handful of steps and reported
    # Success on a 96% error (7 steps here, 60 once the coefficients survive).
    @testset "ExplicitTaylorAdaptiveOrder on Pleiades" begin
        prob = remake(prob_ode_pleiades, tspan = (0.0, 1.0))
        ref = solve(
            prob, ExplicitTaylor(order = Val(6)),
            abstol = 1.0e-12, reltol = 1.0e-12
        )
        sol = solve(
            prob, ExplicitTaylorAdaptiveOrder(min_order = Val(2), max_order = Val(6)),
            abstol = 1.0e-8, reltol = 1.0e-8
        )
        @test SciMLBase.successful_retcode(sol)
        @test length(sol.t) > 30
        @test maximum(abs.(sol.u[end] .- ref.u[end])) < 1.0e-6 * maximum(abs.(ref.u[end]))
    end

    # `get_fsalfirstlast` used to alias the FSAL buffer to `cache.u`, which is
    # `integrator.u`. The auto-dt heuristic writes an extra `f` evaluation into
    # `fsallast`, so the aliasing made that call `f(u, u, p, t)` and corrupted the
    # state before the first step, for any RHS that reads `u` after writing `du`.
    @testset "auto-dt must not write through the state buffer" begin
        function f_alias!(du, u, p, t)
            du[1] = 7.0
            du[2] = 9.0
            s = u[1] + u[2]
            du[1] = -s
            du[2] = s
            return nothing
        end
        u0 = [1.0, 2.0]
        for alg in (ExplicitTaylor(order = Val(4)), ExplicitTaylorAdaptiveOrder())
            integ = init(
                ODEProblem(f_alias!, copy(u0), (0.0, 1.0)), alg,
                abstol = 1.0e-8, reltol = 1.0e-8
            )
            @test !OrdinaryDiffEqTaylorSeries.isfsal(alg)
            @test integ.u == u0
        end
    end

    # End-to-end: with the state corrupted, the auto-dt heuristic saw a NaN and
    # collapsed `dt0` to `dtmin`, so the controller burned hundreds of steps
    # climbing back out (905 before the fix, 587 after, at this tolerance).
    @testset "Auto-dt on Pleiades does not collapse to dtmin" begin
        sol = solve(
            prob_ode_pleiades, ExplicitTaylor(order = Val(6)),
            abstol = 1.0e-10, reltol = 1.0e-10
        )
        @test SciMLBase.successful_retcode(sol)
        @test all(isfinite, sol.u[end])
        @test length(sol.t) < 700
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
            @test length(sol2.t) > 1

            sol8 = solve(
                prob, ExplicitTaylor(order = Val(8)),
                abstol = 1.0e-12, reltol = 1.0e-12
            )
            @test SciMLBase.successful_retcode(sol8)
            @test length(sol8.t) > 1
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
        @test length(sol2.t) > 1

        sol8 = solve(
            prob_oop, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-12, reltol = 1.0e-12
        )
        @test SciMLBase.successful_retcode(sol8)
        @test length(sol8.t) > 1

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
        @test length(sol.t) > 1

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

# Run QA tests (AllocCheck, JET, Aqua) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("qa/allocation_tests.jl")
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end
