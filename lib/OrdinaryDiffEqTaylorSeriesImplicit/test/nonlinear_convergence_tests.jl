using ForwardDiff, LinearAlgebra
import SciMLBase
import OrdinaryDiffEqCore

@testset "Implicit Taylor nonlinear convergence" begin
    ODTI = OrdinaryDiffEqTaylorSeriesImplicit
    square!(du, u, p, t) = (du .= u .^ 2; nothing)
    square_prob = ODEProblem(square!, [1.0], (0.0, 0.5))
    decay!(du, u, p, t) = (du .= -10 .* u; nothing)
    decay_prob = ODEProblem(decay!, [1.0], (0.0, 1.0))

    @testset "Default tolerances and predictor" begin
        for alg in (TaylorGaussA4(), ImplicitTaylor1A2())
            i = init(square_prob, alg; dt = 0.05)
            @test i.cache.κ == 0.01
            @test alg.extrapolant === :extrapolate
        end
        custom = init(square_prob, TaylorGaussA4(κ = 1e-5); dt = 0.05)
        @test custom.cache.κ == 1e-5
    end

    @testset "Stage solutions match the Padé formulas" begin
        problem = decay_prob
        # [2/2] (Taylor-Gauss A4) and [1/2] (Taylor-Radau LA3) on y' = -10y.
        formulas = (
            (TaylorGaussA4, z -> (1 + z / 2 + z^2 / 12) / (1 - z / 2 + z^2 / 12)),
            (TaylorRadauLA3, z -> (1 + z / 3) / (1 - 2z / 3 + z^2 / 6)),
        )
        for (make, pade) in formulas, extrapolant in (:constant, :linear, :polynomial)
            for adaptive in (false, true), h in (0.01, 1.0)
                i = init(problem, make(; extrapolant); dt = h, adaptive)
                out = similar(i.u)
                success = ODTI.implicit_taylor_stage!(i, i.cache, 0.0, h, [1.0], out)
                @test success
                @test only(out) ≈ pade(-10h)
                @test i.cache.status == OrdinaryDiffEqCore.Convergence
            end
        end
    end

    @testset "Eta still prevents premature acceptance in slow convergence" begin
        # Backward Euler for u' = u^2, u_n = 1, h = 0.24 has the nearby
        # root 5/3. With the Jacobian frozen at u_n (constant predictor) the
        # contraction tends to 8/13. At iteration 6 (κ = 0.01) or 20
        # (κ = 1e-5), the correction is below κ, but its remaining tail is
        # still above κ.
        for (κ, early_stop) in ((0.01, 6), (1e-5, 20)), maxiters in (early_stop, 40)
            i = init(
                square_prob, ImplicitTaylor1L1(; κ, maxiters, extrapolant = :constant);
                dt = 0.24, adaptive = false, abstol = 1.0, reltol = 0.0
            )
            out = similar(i.u)
            success = ODTI.implicit_taylor_stage!(
                i, i.cache, 0.0, 0.24, [1.0], out
            )
            if maxiters == early_stop
                @test !success
                @test i.force_stepfail
            else
                @test success
                @test i.cache.ηold > 1
                @test i.cache.iter > early_stop
                @test abs(only(out) - 5 / 3) < i.cache.κ
            end
        end
    end

    @testset "Equal correction norms do not certify convergence" begin
        # Backward Euler for u' = 1 - u^2 from u_n = 0, h = 1 has
        # F(u) = u^2 + u - 1. Freezing F'(0) = 1 makes the iterates
        # alternate between 0 and 1. Equal correction norms are a cycle,
        # not floating-point stagnation or convergence to the nearby root.
        cycling!(du, u, p, t) = (du .= 1 .- u .^ 2; nothing)
        problem = ODEProblem(cycling!, [0.0], (0.0, 1.0))
        for adaptive in (false, true)
            i = init(problem, ImplicitTaylor1L1(extrapolant = :constant, maxiters = 3);
                dt = 1.0, adaptive, abstol = 1.0, reltol = 0.0)
            @test !ODTI.implicit_taylor_stage!(i, i.cache, 0.0, 1.0, [0.0], i.u)
            @test i.force_stepfail
        end
    end

    @testset "Nonlinear failure refreshes the retry linearization" begin
        # On u' = u^2, c_k(u) = u^(k+1), so J_1 = 2u. The larger
        # step exhausts a short Newton budget. Reducing h changes the
        # linear predictor, and the retry must evaluate J_k there even
        # though the start time and state are unchanged.
        for make in (TaylorGaussA4, TaylorRadauLA3, TaylorRadauLA5,
                TaylorLobattoLA4, ImplicitTaylor1L1, ImplicitTaylor1A2,
                ImplicitTaylor2A4)
            alg = make(extrapolant = :linear, maxiters = 3)
            options = (abstol = 1e-6, reltol = 1e-6)
            i = init(square_prob, alg; dt = 0.4, options...)
            c = i.cache
            out = similar(i.u)
            @test !ODTI.implicit_taylor_stage!(i, c, 0.0, 0.4, [1.0], out)
            @test i.force_stepfail
            @test !c.jac_valid

            h = 0.01
            njacs = i.stats.njacs
            i.force_stepfail = false
            @test ODTI.implicit_taylor_stage!(i, c, 0.0, h, [1.0], out)
            # The retry refreshes its predictor Jacobian, and convergence
            # validation also evaluates it at the final implicit point.
            @test i.stats.njacs >= njacs + 2
            endpoint_jacobian = size(c.Jk, 3) == 1 ? c.Jk[1, 1, 1] : c.Mtmp[1, 1]
            @test endpoint_jacobian ≈ 2 * only(c.uintermediate)

            # The failed larger step must not contaminate a fresh solve
            # at the reduced step size, including complex μ-Taylor stages.
            fresh = init(square_prob, alg; dt = h, options...)
            expected = similar(out)
            @test ODTI.implicit_taylor_stage!(
                fresh, fresh.cache, 0.0, h, [1.0], expected
            )
            @test out ≈ expected rtol = 10eps(Float64)
        end
    end

    @testset "Error rejection preserves coefficient Jacobian reuse" begin
        # This linear stage converges, but its error estimate rejects the
        # large step. Retrying it must not be treated as a Newton failure.
        i = init(decay_prob, TaylorRadauLA3(); dt = 0.5, abstol = 1e-8, reltol = 1e-6)
        c = i.cache
        out = similar(i.u)
        @test ODTI.implicit_taylor_step!(i, c, 0.0, 0.5, [1.0], out)
        @test OrdinaryDiffEqCore.get_EEst(i) > 1
        @test c.jac_valid && !i.force_stepfail
        njacs = i.stats.njacs
        h = 0.1
        @test ODTI.implicit_taylor_step!(i, c, 0.0, h, [1.0], out)
        # The retry reuses the starting Jacobian; the one extra call is
        # the mandatory check at the candidate stage.
        @test i.stats.njacs == njacs + 1
        z = -10h
        @test only(out) ≈ (1 + z / 3) / (1 - 2z / 3 + z^2 / 6)
    end

    @testset "Adaptive mu-Taylor accuracy" begin
        for make in (ImplicitTaylor1A2, ImplicitTaylor2A4)
            errors = Float64[]
            for tol in (1e-3, 1e-6)
                sol = solve(square_prob, make(); abstol = tol / 100, reltol = tol)
                @test SciMLBase.successful_retcode(sol)
                err = maximum(abs(only(u) - 1 / (1 - t)) for (t, u) in zip(sol.t, sol.u))
                @test err < 2tol
                push!(errors, err)
            end
            @test errors[2] < errors[1] / 100
        end
    end

    @testset "Small initial contraction ratios on HIRES" begin
        hires = ODEProblemLibrary.prob_ode_hires
        rhs! = SciMLBase.unwrapped_f(hires.f)
        # Evaluate the [2/2] and [1/2] equations independently of the
        # generated Taylor polynomials and the solver's stopping rule.
        function reference_stage(uprev, t, h, guess, embedded)
            function rate(v)
                du = similar(v)
                rhs!(du, v, hires.p, t)
                return du
            end
            second(v) = ForwardDiff.jacobian(rate, v) * rate(v)
            a1, a2, b1, b2 = embedded ?
                (1 / 3, 0.0, -2 / 3, 1 / 6) : (1 / 2, 1 / 12, -1 / 2, 1 / 12)
            rhs = a1 * h * rate(uprev) + a2 * h^2 * second(uprev)
            residual(z) = z + b1 * h * rate(uprev + z) +
                b2 * h^2 * second(uprev + z) - rhs
            z = guess - uprev
            for _ in 1:30
                correction = ForwardDiff.jacobian(residual, z) \ residual(z)
                z -= correction
                if norm(correction, Inf) < 1e-14
                    return uprev + z
                end
            end
            error("Independent HIRES stage solve did not converge")
        end
        # Recorded starts of accepted steps that previously underestimated
        # the remaining error in both the main and embedded solves.
        cases = (
            (
                t = 108.06757323055373, h = 4.702933562660071, κ = 1e-5,
                u = [
                    0.004376219796853137, 0.0008556386996123897,
                    0.0007689288265462751, 0.007558171943033635,
                    0.1273954175820801, 0.5099217674805134,
                    0.0056286476291058435, 7.13523708941554e-5,
                ],
            ),
            (
                t = 218.10456007814102, h = 1.2656587632349527, κ = 0.01,
                u = [
                    0.0024168654663418625, 0.00047271809059855075,
                    0.0003860519682369947, 0.004128578038121805,
                    0.0598105764371007, 0.2357590001276736,
                    0.005547752763925982, 0.0001522472360740147,
                ],
            ),
        )
        # The [1/2] equation (embedded in Taylor-Gauss A4) is solved as the
        # main stage of Taylor-Radau LA3. Fixed steps also exercise the
        # stopping rule.
        for case in cases, embedded in (false, true), adaptive in (false, true)
            make = embedded ? TaylorRadauLA3 : TaylorGaussA4
            i = init(
                hires, make(κ = case.κ, extrapolant = :polynomial);
                dt = case.h, abstol = 1e-6, reltol = 1e-3, adaptive
            )
            out = similar(case.u)
            success = ODTI.implicit_taylor_stage!(
                i, i.cache, case.t, case.h, case.u, out
            )
            @test success
            if success
                reference = reference_stage(case.u, case.t, case.h, out, embedded)
                scale = 1e-6 .+ 1e-3 .* max.(abs.(case.u), abs.(out))
                actual_tail = norm((out - reference) ./ scale) / sqrt(length(out))
                @test actual_tail < case.κ
                @test i.cache.iter > 2
                @test i.cache.ηold > 0
            end
        end
    end

    # Shared by the testsets below: every new (function, method) pair costs a
    # symbolic code generation, so they reuse this problem where they can.
    stiff(u, p, t) = -1000 * (u - cos(t)) - sin(t)
    stiff_prob = ODEProblem(stiff, 1.0, (0.0, 1.0))

    @testset "Adaptive stiff nonautonomous scalar problem" begin
        problem = stiff_prob
        for alg in (TaylorGaussA4(), TaylorRadauLA3())
            errors = Float64[]
            for tol in (1e-3, 1e-7)
                sol = solve(problem, alg; dt = 0.1, abstol = tol / 100, reltol = tol)
                @test SciMLBase.successful_retcode(sol)
                push!(errors, maximum(abs.(sol.u .- cos.(sol.t))))
            end
            @test errors[2] < errors[1]
            @test errors[2] < 1e-6
        end
    end

    @testset "Step size controllers" begin
        @test TaylorGaussA4().controller === :Predictive
        @test_throws ArgumentError TaylorGaussA4(controller = :I)
        for (controller, cachetype) in (
                (:Predictive, OrdinaryDiffEqCore.PredictiveControllerCache),
                (:PI, OrdinaryDiffEqCore.PIControllerCache),
            ), make in (TaylorGaussA4, TaylorRadauLA3)
            i = init(stiff_prob, make(; controller); dt = 0.1, abstol = 1e-9, reltol = 1e-7)
            @test i.controller_cache isa cachetype
            sol = solve!(i)
            @test SciMLBase.successful_retcode(sol)
            @test maximum(abs.(sol.u .- cos.(sol.t))) < 1e-5
        end
    end

    @testset "Second-order Padé methods estimate with third divided differences" begin
        # u' = -2(u - cos t) - sin t has the solution cos t for u(0) = 1.
        smooth(u, p, t) = -2 * (u - cos(t)) - sin(t)
        problem = ODEProblem(smooth, 1.0, (0.0, 20.0))
        for alg in (TaylorGaussA2(), TaylorLobattoLA2())
            @test ODTI.uses_divided3(alg)
            @test !ODTI.uses_embedded(alg)
            @test ODTI.alg_adaptive_order(alg) == 2
            steps = Int[]
            errors = Float64[]
            for tol in (1e-3, 1e-6)
                sol = solve(problem, alg; abstol = tol, reltol = tol)
                @test SciMLBase.successful_retcode(sol)
                push!(steps, sol.stats.naccept)
                push!(errors, maximum(abs.(sol.u .- cos.(sol.t))))
            end
            # Errors per step accumulate over the 20 time units (by about
            # 1/(2h) on this problem); the global error still scales with
            # the tolerance.
            @test errors[2] < 100 * 1e-6
            @test errors[2] / errors[1] < 1e-2
            # An estimate of the right order grows the step count like
            # tol^(-1/3): 1000× tighter means about 10× more steps, not the
            # ~30× of a first-order estimate.
            @test 4 < steps[2] / steps[1] < 16
        end
        # The embedded estimate is still used from order three on.
        @test ODTI.uses_embedded(TaylorRadauLA3())
        @test ODTI.alg_adaptive_order(TaylorRadauLA3()) == 2
    end

    @testset "Residual estimates are filtered through (I - hJ/(P+Q))^Q" begin
        # On y' = λy the embedded residual at the computed state is
        # r1 = (D1 N - N1 D)/D · y0, with N/D and N1/D1 the Padé
        # approximants of the method and of its embedded formula, and the
        # estimate divides it by the filter (1 - hλ/(P+Q))^Q rather than by
        # the Newton matrix D(hλ) (see implicit_taylor_error!).
        pade_poly(P, Q) = (
            [factorial(P) * factorial(P + Q - k) / (factorial(P + Q) * factorial(k) * factorial(P - k)) for k in 0:P],
            [(-1)^k * factorial(Q) * factorial(P + Q - k) / (factorial(P + Q) * factorial(k) * factorial(Q - k)) for k in 0:Q],
        )
        polyval(c, w) = sum(c[k + 1] * w^k for k in 0:(length(c) - 1))
        methods = (
            (TaylorGaussA4(), 2, 2), (TaylorRadauLA3(), 1, 2), (TaylorLobattoLA4(), 1, 3),
            (ImplicitTaylor(order = Val(2), order_q = Val(1)), 2, 1),
        )
        for (alg, P, Q) in methods, h in (0.01, 0.3, 3.0)
            P1, Q1 = ODTI.embedded_degrees(P, Q)
            N, D = pade_poly(P, Q)
            N1, D1 = pade_poly(P1, Q1)
            w = -10h
            r1 = (polyval(D1, w) * polyval(N, w) - polyval(N1, w) * polyval(D, w)) / polyval(D, w)
            expected = r1 / (1 - w / (P + Q))^Q
            i = init(decay_prob, alg; dt = h, abstol = 1e-8, reltol = 1e-6)
            out = similar(i.u)
            factorizations = i.stats.nw
            @test ODTI.implicit_taylor_step!(i, i.cache, 0.0, h, [1.0], out)
            u = only(out)
            @test u ≈ polyval(N, w) / polyval(D, w)
            scale = 1e-8 + 1e-6 * max(1.0, abs(u))
            @test OrdinaryDiffEqCore.get_EEst(i) ≈ abs(expected) / scale rtol = 1e-6
            # The Newton matrix, and the filter unless it is the Newton matrix.
            @test i.stats.nw - factorizations == (Q == 1 ? 1 : 2)
        end
    end

    @testset "Small dense filters preserve repeated solves and solver overrides" begin
        LS = ODTI.LinearSolve
        small = init(ODEProblem(decay!, ones(8), (0.0, 1.0)), TaylorRadauLA5()).cache
        expected = LS.defaultalg(small.W, vec(small.tmp), LS.OperatorAssumptions(true))
        @test small.linsolve_filter.b isa AbstractMatrix
        @test small.linsolve_filter.alg.alg == expected.alg
        for T in (Float32, Float64, ComplexF64),
                solver in (nothing, LS.LUFactorization(), LS.SimpleLUFactorization())
            problem = ODEProblem(decay!, ones(T, 20), (0.0, 1.0))
            i = init(problem, TaylorRadauLA5(linsolve = solver); dt = 0.1)
            c = i.cache
            solver === nothing || @test c.linsolve_filter.alg === solver
            # A coupled, stiff Jacobian exercises the filter independently of
            # the Newton iteration and its choice of factorization backend.
            J = zeros(T, 20, 20)
            for j in 1:20
                J[j, j] = -10j
                j < 20 && (J[j, j + 1] = 2)
            end
            c.Jk[:, :, 1] .= J
            rhs = T.(1:20)
            for (h, fresh) in ((0.1, true), (0.1, false), (0.3, true))
                W = Matrix{T}(I, 20, 20) - T(h / 5) * J
                reference = copy(rhs)
                for _ in 1:3
                    reference = W \ reference
                end
                c.tmp .= rhs
                nw, ns = i.stats.nw, i.stats.nsolve
                ODTI.filter_error!(i, c, c.tmp, 3, -3 / 5, h)
                @test c.tmp ≈ reference rtol = 100eps(real(T))
                @test i.stats.nw - nw == fresh
                @test i.stats.nsolve - ns == 3
            end
            if c.linsolve_filter.b isa AbstractMatrix
                # The two columns have independent right-hand sides and
                # share the same factorization and repeated filter solves.
                h = 0.1
                W = Matrix{T}(I, 20, 20) - T(h / 5) * J
                first, second = copy(rhs), reverse(rhs)
                expected_first, expected_second = copy(first), copy(second)
                for _ in 1:3
                    expected_first = W \ expected_first
                    expected_second = W \ expected_second
                end
                ns = i.stats.nsolve
                ODTI.filter_error!(i, c, first, 3, -3 / 5, h, nothing, second)
                @test first ≈ expected_first rtol = 100eps(real(T))
                @test second ≈ expected_second rtol = 100eps(real(T))
                @test i.stats.nsolve - ns == 6
            end
        end
    end

    @testset "Linear solvers that return the solution in their own buffer" begin
        # SimpleLUFactorization leaves the solution in its workspace rather
        # than in the cache's `u`, where the step reads it; the cache detects
        # this at initialization and copies.
        LS = ODTI.LinearSolve
        for make in (TaylorGaussA4, TaylorRadauLA5)
            i = init(stiff_prob, make(linsolve = LS.SimpleLUFactorization()); dt = 0.1, abstol = 1e-9, reltol = 1e-7)
            c = i.cache isa ODTI.ImplicitTaylorConstantCache ? i.cache.stage_cache : i.cache
            @test c.detached_solution
            sol = solve!(i)
            @test SciMLBase.successful_retcode(sol)
            @test maximum(abs.(sol.u .- cos.(sol.t))) < 1e-5
        end
        c = init(stiff_prob, TaylorGaussA4(); dt = 0.1).cache.stage_cache
        @test !c.detached_solution && !c.detached_filter
    end

    @testset "Factorizations are kept while the matrices change little" begin
        # y' = -10y with Lobatto LA4: the Newton matrix D(hJ) is kept while
        # it changes by less than 1/20, the filter factor I - hJ/4 while it
        # changes by less than 1/5 of max(1, its norm). The estimate then
        # uses the kept factor, which the closed form shows.
        pade_poly(P, Q) = (
            [factorial(P) * factorial(P + Q - k) / (factorial(P + Q) * factorial(k) * factorial(P - k)) for k in 0:P],
            [(-1)^k * factorial(Q) * factorial(P + Q - k) / (factorial(P + Q) * factorial(k) * factorial(Q - k)) for k in 0:Q],
        )
        polyval(c, w) = sum(c[k + 1] * w^k for k in 0:(length(c) - 1))
        N, D = pade_poly(1, 3)
        N1, D1 = pade_poly(1, 2)
        estimate(w, w_filter, y) = abs(
            (polyval(D1, w) * polyval(N, w) - polyval(N1, w) * polyval(D, w)) / polyval(D, w) /
                (1 - w_filter / 4)^3 * y
        )
        i = init(decay_prob, TaylorLobattoLA4(); dt = 1.0, abstol = 1e-8, reltol = 1e-6)
        c = i.cache
        t = 0.0
        y = 1.0
        h_filter = 0.0
        # (h, Newton matrix refactorized, filter refactorized): an unchanged
        # h changes neither; h = 1.1 changes D(hJ) by 26% but the filter
        # factor 2.5 of the first step by 0.25 < 0.5; h = 1.4 changes both.
        for (h, newton, filter) in ((1.0, true, true), (1.0, false, false), (1.1, true, false),
                (1.4, true, true), (1.4, false, false))
            out = similar(i.u)
            factorizations = i.stats.nw
            @test ODTI.implicit_taylor_step!(i, c, t, h, [y], out)
            @test i.stats.nw - factorizations == newton + filter
            filter && (h_filter = h)
            u = only(out)
            scale = 1e-8 + 1e-6 * max(abs(y), abs(u))
            @test OrdinaryDiffEqCore.get_EEst(i) ≈ estimate(-10h, -10h_filter, y) / scale rtol = 1e-6
            t += h
            y = u
        end
    end

    @testset "Degree-3 methods reject roots on the wrong part of the slow manifold" begin
        # With the (0, 3) embedded formula, Lobatto LA4 at reltol 1e-3 took a
        # step of 94 time units over HIRES's slow phase, accepted with
        # EEst = 0.003 and an error of 10 tolerances: the step equation's
        # roots far beyond the nonlinearity limit lie on the slow manifold
        # and a formula of the same implicit degree shares them. The
        # (1, 2) embedded formula does not; its residual at such a root
        # (185 tolerances for a step of 112 time units) is then let
        # through by the filter of the estimate, where the Newton matrix
        # had damped it by 10⁴ (implicit_taylor_error!).
        @test ODTI.embedded_degrees(1, 3) == (1, 2)
        @test ODTI.embedded_degrees(2, 3) == (2, 2)
        @test ODTI.embedded_degrees(1, 2) == (0, 2)
        @test ODTI.embedded_degrees(0, 2) == (0, 1)
        hires = ODEProblemLibrary.prob_ode_hires
        reference = solve(hires, Rodas5P(); abstol = 1e-11, reltol = 1e-11, dense = true)
        scaled_error(u, e) = sqrt(sum(abs2, (u .- e) ./ (1e-6 .+ 1e-3 .* abs.(e))) / length(u))
        for alg in (TaylorLobattoLA4(), TaylorRadauLA5())
            sol = solve(hires, alg; abstol = 1e-6, reltol = 1e-3)
            @test SciMLBase.successful_retcode(sol)
            worst = maximum(scaled_error(u, reference(t)) for (t, u) in zip(sol.t, sol.u))
            @test worst < 5
        end
        # Single Lobatto LA4 steps from exact states over the slow phase,
        # with step sizes the Newton-matrix-filtered estimate accepted. From
        # most of these starts Newton converges to roots on the slow
        # manifold with errors of 1.5-60 tolerances, which the Newton
        # matrix damped to EEst < 0.03; near the nonlinearity limit the
        # iteration is sensitive to its start, so which root it reaches, or
        # whether it fails (which rejects the step as well), is not
        # asserted. Whatever root it reaches, the estimate must not be blind
        # to its error (a far root may be overestimated: the residual only
        # measures the error to first order) and must be on its side of 1,
        # and at least two starts must be rejected, either by the nonlinear
        # solve or by the error estimate. The estimator is also exercised
        # independently at a spurious root below.
        # The last start is a legitimate step (error 0.2 tolerances) that
        # must be accepted. The extrapolating predictor gets two earlier
        # exact states, with the growth of the old run.
        cases = (
            (98.0, 112.0, :polynomial), (50.0, 25.0, :polynomial), (98.0, 55.0, :polynomial),
            (150.0, 25.0, :polynomial), (150.0, 55.0, :polynomial), (210.0, 55.0, :extrapolate),
            (98.0, 55.0, :extrapolate),
        )
        rejected_roots = 0
        for (t0, h, extrapolant) in cases
            u0 = reference(t0)
            problem = SciMLBase.remake(hires; u0 = copy(u0), tspan = (t0, t0 + h))
            i = init(problem, TaylorLobattoLA4(; extrapolant, maxiters = 50); dt = h, abstol = 1e-6, reltol = 1e-3)
            c = i.cache
            if extrapolant === :extrapolate
                h1, h2 = h / 2, h / 2.4
                i.tprev = t0 - h1
                c.tprev2 = t0 - h1 - h2
                copyto!(c.uprev2, reference(t0 - h1))
                copyto!(c.uprev3, reference(t0 - h1 - h2))
                i.success_iter = 2
            end
            out = similar(u0)
            converged = ODTI.implicit_taylor_step!(i, c, t0, h, copy(u0), out)
            legitimate = (t0, h, extrapolant) == (98.0, 55.0, :extrapolate)
            legitimate && @test converged
            if !converged
                @test i.force_stepfail
                rejected_roots += 1
                continue
            end
            err = scaled_error(out, reference(t0 + h))
            EEst = OrdinaryDiffEqCore.get_EEst(i)
            @test EEst > err / 10
            @test (EEst > 1) == (err > 1)
            legitimate && @test err < 1
            rejected_roots += EEst > 1
        end
        @test rejected_roots >= 2

        # A spurious root of the scalar [1/3] stage equation. Test the
        # estimator directly, independently of whether the nonlinear solver
        # now refuses to approach this root from its usual predictor.
        h = 0.1
        a, b = ODTI.pade_increment_coefficients(Val(1), Val(3))
        root = big"30"
        for _ in 1:40
            residual = root - 1 - BigFloat(a[2]) * BigFloat(h)
            derivative = one(root)
            for k in 1:3
                weight = BigFloat(b[k + 1]) * BigFloat(h)^k
                residual += weight * root^(k + 1)
                derivative += (k + 1) * weight * root^k
            end
            correction = residual / derivative
            root -= correction
            abs(correction) < big"1e-40" && break
        end
        @test root > 10
        @test abs(root - 1 - BigFloat(a[2]) * BigFloat(h) +
            sum(BigFloat(b[k + 1]) * BigFloat(h)^k * root^(k + 1) for k in 1:3)) < big"1e-35"
        i = init(square_prob, TaylorLobattoLA4(); dt = h, abstol = 1e-6, reltol = 1e-3)
        c = i.cache
        c.polynomial_prev(c.rhs, c.rhs_embedded, c.predictor, [1.0], 0.0, h)
        c.uintermediate[1] = Float64(root)
        c.utilde[1] = Float64(root) - 1
        c.polynomial_jac(c.tmp, c.tmp_embedded, c.Jk, c.uintermediate, h, h)
        c.Mtmp .= @view c.Jk[:, :, 1]
        ODTI.implicit_taylor_error!(i, c, 0.0, h, [1.0], c.uintermediate)
        @test OrdinaryDiffEqCore.get_EEst(i) > 1
    end

    @testset "The step size cap is probed with backoff" begin
        c = init(decay_prob, TaylorLobattoLA4(); dt = 1.0).cache
        c.dt_cap = 1.0
        # An accepted step with a fast iteration, as the controller reports it.
        fast!(dt) = ODTI.newton_feedback!(c, 1.0, dt)
        # One fast step at the cap raises it at once.
        fast!(1.0)
        @test c.dt_cap ≈ 1.25 && c.cap_raised
        # The step at the raised cap fails: the cap returns to 4/5 of the
        # failed size and the next raise waits for two fast steps at it.
        ODTI.record_newton_failure!(c, 1.25)
        @test c.dt_cap ≈ 1.0 && c.cap_probe == 2 && !c.cap_raised
        c.newton_failed = false
        fast!(0.5)
        fast!(1.0)
        @test c.dt_cap ≈ 1.0 && c.cap_streak == 1
        fast!(1.0)
        @test c.dt_cap ≈ 1.25 && c.cap_raised
        # Another failure there: four fast steps before the next probe.
        ODTI.record_newton_failure!(c, 1.25)
        c.newton_failed = false
        @test c.cap_probe == 4
        for _ in 1:3
            fast!(1.0)
        end
        @test c.dt_cap ≈ 1.0
        fast!(1.0)
        @test c.dt_cap ≈ 1.25 && c.cap_raised
        # A fast step at the raised cap: the limit receded, probe every step again.
        fast!(1.25)
        @test c.cap_probe == 1 && c.dt_cap ≈ 1.25^2
        # A failure below the cap is a new limit: probe every step.
        c.cap_probe = 4
        ODTI.record_newton_failure!(c, 0.5)
        @test c.cap_probe == 1 && c.dt_cap ≈ 0.4
        # Neither a failed nor a slow step counts towards a raise.
        c.dt_cap = 1.0
        c.newton_failed = true
        fast!(1.0)
        @test c.dt_cap ≈ 1.0 && c.cap_streak == 0
        c.newton_limit = true
        fast!(1.0)
        @test c.dt_cap ≈ 1.0 && c.cap_streak == 0
    end

    @testset "Extrapolating predictor reproduces quadratics" begin
        poly(t) = [1 + t + 3t^2]
        for (alg, pade) in ((TaylorGaussA4(), true), (ImplicitTaylor1A2(), false))
            i = init(square_prob, alg; dt = 0.05)
            c = i.cache
            t, tprev, tprev2 = 0.7, 0.4, 0.3
            c.tprev2 = tprev2
            i.tprev = tprev
            copyto!(c.uprev3, poly(tprev2))
            copyto!(c.uprev2, poly(tprev))
            uprev = poly(t)
            h = 0.5
            distance = pade ? h : alg.μ * h
            # Two earlier states: linear extrapolation.
            i.success_iter = 1
            ODTI.taylor_initial_guess!(i, c, alg, t, h, uprev, pade)
            @test c.utilde ≈ (poly(t) - poly(tprev)) * (distance / (t - tprev))
            # Three earlier states: the quadratic is reproduced exactly.
            i.success_iter = 2
            ODTI.taylor_initial_guess!(i, c, alg, t, h, uprev, pade)
            @test c.utilde ≈ poly(t + distance) - uprev
            # A step more than 3× the previous one falls back to the linear guess.
            ODTI.taylor_initial_guess!(i, c, alg, t, 2.0, uprev, pade)
            @test c.utilde ≈ (poly(t) - poly(tprev)) * ((pade ? 2.0 : alg.μ * 2.0) / (t - tprev))
        end
    end
end

include("heldout_tests.jl")
