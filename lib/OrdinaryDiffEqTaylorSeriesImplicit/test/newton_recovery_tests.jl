@testset "Implicit Taylor Newton recovery" begin
    ODTI = OrdinaryDiffEqTaylorSeriesImplicit
    square!(du, u, p, t) = (du .= u .^ 2; nothing)
    problem = ODEProblem(square!, [1.0], (0.0, 0.5))

    @testset "A stale LU retry keeps the current predictor Jacobians" begin
        # M = [1 100; 0 1], B = [1 100; .02 1]. Their row change is
        # below 5%, but I-B\M has eigenvalue 2, so a retained B diverges.
        # The linear equation has an exact backward-Euler solution and its
        # predictor derivatives remain valid when the solve restarts there.
        A = [0.0 -1000.0; 0.0 0.0]
        linear!(du, u, p, t) = (mul!(du, A, u); nothing)
        for scale in (1.0e-6, 1.0, 1.0e6), controller in (:PI, :Predictive),
                prior_failure in (false, true)
            u0 = scale .* [1.0, 0.01]
            p = ODEProblem(linear!, u0, (0.0, 1.0))
            i = init(
                p, ImplicitTaylor1L1(; extrapolant = :constant, controller);
                dt = 0.1, abstol = 1.0e-8scale, reltol = 0.0
            )
            c = i.cache
            kept = [1.0 100.0; 0.02 1.0]
            copyto!(c.J, kept)
            copyto!(c.Mkept, kept)
            ODTI.matrix_row_norms!(c.Mrow_norm, kept)
            # Install the kept factorization the way the solver stores it.
            if c.direct_lu
                ODTI.dense_lu!(c.J, c.pivJ)
            else
                ODTI.dolinsolve(i, c.linsolve; A = c.J, b = zeros(2), linu = c.linsolve.u)
            end
            c.lu_valid = true
            # Force the first attempt to evaluate the current Jk. Its retry
            # must keep that result rather than relying on the constant flag.
            c.jac_valid = c.constant_jacobian = false
            c.newton_failed = prior_failure
            starts = Vector{Float64}[]
            original = c.polynomial_jac
            c.polynomial_jac = typeof(original)(
                (b, b1, Jk, u, t, h) -> begin
                    push!(starts, copy(u))
                    original(b, b1, Jk, u, t, h)
                    nothing
                end
            )
            @test ODTI.implicit_taylor_stage!(i, c, 0.0, 0.1, u0, i.u)
            @test c.nrestart == 1
            @test c.ncoefjac == 1
            @test length(starts) == 2 # predictor, then the independent endpoint check
            @test first(starts) == u0
            @test last(starts) ≈ i.u
            @test c.nendpoint_refresh == c.nvalidation_fail == 0
            @test c.ncheckjac == 1
            @test c.ncoefjac + c.ncheckjac == i.stats.njacs
            @test norm(i.u - scale .* [0.0, 0.01], Inf) < c.κ * 1.0e-8scale
            # Two discarded chord solves, then at most two fresh corrections;
            # no third old-LU correction or auxiliary probe solve is needed.
            @test c.iter == i.stats.nnonliniter <= 4
            @test i.stats.nsolve <= 5 # four corrections + degree-one check
            @test !c.recovery_since_accept && !c.newton_limit
            @test c.ηold == 1 # No rate was measured in the successful two-correction attempt.
            @test c.newton_failed == prior_failure && iszero(c.dt_cap)
            OrdinaryDiffEqCore.set_EEst!(i, 0.01)
            next_h = OrdinaryDiffEqCore.step_accept_controller!(i, i.alg, 0.2)
            @test prior_failure ? next_h <= 0.1 : next_h > 0.1
            @test !c.newton_failed && !c.recovery_since_accept
        end
    end

    @testset "Stale-matrix failure restarts at the original predictor" begin
        i = init(
            problem, ImplicitTaylor1L1(extrapolant = :constant, maxiters = 40);
            dt = 0.1, abstol = 1.0e-6, reltol = 0.0
        )
        c = i.cache
        # Simulate a retry retaining the derivative from u=-10. Reassembling
        # and factorizing it at this h does not make the linearization fresh.
        c.Jk .= -20
        c.jac_valid = true
        starts = Float64[]
        original = c.polynomial_jac
        c.polynomial_jac = typeof(original)(
            (b, b1, Jk, u, t, h) -> begin
                push!(starts, only(u))
                original(b, b1, Jk, u, t, h)
                nothing
            end
        )
        @test ODTI.implicit_taylor_stage!(i, c, 0.0, 0.1, [1.0], i.u)
        @test c.nrestart == 1
        @test first(starts) == 1.0
        physical_root = 2 / (1 + sqrt(1 - 4 * 0.1))
        @test abs(only(i.u) - physical_root) < c.κ * 1.0e-6
        @test c.iter == i.stats.nnonliniter
        @test c.iter <= 2i.alg.maxiters
        @test c.ncoefjac + c.ncheckjac == i.stats.njacs
        @test !c.recovery_since_accept && !c.newton_limit
        @test !c.newton_failed && iszero(c.dt_cap)
    end

    @testset "A slow fresh retry retains feedback through an error retry" begin
        # F(z)=.2z-.2z²-.018 at h=.2, with physical root z=.1. The
        # fresh predictor matrix is .2, the endpoint derivative is .16,
        # and the chord contraction is .2 (eta=.25), above the fast cutoff.
        shifted_square!(du, u, p, t) = (du .= p .* (4 .* (u ./ p .- 1) .+ (u ./ p .- 1) .^ 2 .+ 0.09); nothing)
        for scale in (1.0e-6, 1.0, 1.0e6)
            p = ODEProblem(shifted_square!, [scale], (0.0, 1.0), scale)
            i = init(
                p, ImplicitTaylor1L1(extrapolant = :constant, maxiters = 40);
                dt = 0.2, abstol = 1.0e-9scale, reltol = 0.0
            )
            c = i.cache
            c.Jk .= 20 # Borrowed derivatives must be reevaluated at the predictor.
            c.jac_valid = true
            @test ODTI.implicit_taylor_stage!(i, c, 0.0, 0.2, [scale], i.u)
            @test c.nrestart == 1 && c.ncoefjac == 1
            @test c.nendpoint_refresh == c.nvalidation_fail == 0
            @test c.ηold > ODTI.fast_convergence_cutoff(i.alg)
            @test c.recovery_since_accept && c.newton_limit
            @test abs(only(i.u) - 1.1scale) < c.κ * 1.0e-9scale
            # A smaller error retry must not erase the slow fresh attempt.
            @test ODTI.implicit_taylor_stage!(i, c, 0.0, 0.02, [scale], i.u)
            @test c.recovery_since_accept && c.newton_limit
            @test ODTI.newton_feedback!(c, 0.2, 0.02) == 1
            @test !c.recovery_since_accept && !c.newton_failed
        end
    end

    @testset "Endpoint recovery is bounded and retains its step size feedback" begin
        # Backward Euler on u'=u^2 has two roots for h<1/4. With a
        # constant predictor the matrix changes within the step, triggering
        # validation and one fresh endpoint correction. Very loose scaling
        # nominates a candidate earlier, before the refreshed matrix has
        # stabilized; a second failed check must reject, not polish again.
        for (h, abstol, expected) in ((0.2, 1.0, true), (0.24, 100.0, false))
            i = init(
                problem, ImplicitTaylor1L1(extrapolant = :constant, maxiters = 40);
                dt = h, abstol, reltol = 0.0
            )
            c = i.cache
            @test ODTI.implicit_taylor_stage!(i, c, 0.0, h, [1.0], i.u) == expected
            @test c.nendpoint_refresh == 1
            @test c.nrestart == 0
            @test c.ncoefjac == 1 && c.ncheckjac == 2
            @test c.nvalidation_fail == (expected ? 1 : 2)
            @test c.iter == i.stats.nnonliniter <= i.alg.maxiters
            @test c.recovery_since_accept
            if expected
                @test c.newton_limit
                @test !c.newton_failed && iszero(c.dt_cap)
                physical_root = 2 / (1 + sqrt(1 - 4h))
                @test abs(only(i.u) - physical_root) < c.κ * abstol
            else
                @test i.force_stepfail && c.newton_failed
            end
        end
        i = init(
            problem, ImplicitTaylor1L1(extrapolant = :constant, maxiters = 3);
            dt = 0.2, abstol = 1.0, reltol = 0.0
        )
        @test !ODTI.implicit_taylor_stage!(i, i.cache, 0.0, 0.2, [1.0], i.u)
        @test i.cache.iter <= i.alg.maxiters
        @test i.cache.nendpoint_refresh == 0
    end

    @testset "Fresh-matrix divergence is not polished at the bad iterate" begin
        cycle!(du, u, p, t) = (du .= 1 .- u .^ 2; nothing)
        p = ODEProblem(cycle!, [0.0], (0.0, 1.0))
        i = init(
            p, ImplicitTaylor1L1(extrapolant = :constant, maxiters = 10);
            dt = 1.0, abstol = 1.0, reltol = 0.0
        )
        c = i.cache
        @test !ODTI.implicit_taylor_stage!(i, c, 0.0, 1.0, [0.0], i.u)
        @test c.nrestart == c.nendpoint_refresh == c.ncheckjac == 0
        @test c.ncoefjac == 1
        @test c.iter == 3
        @test c.nfilter_recovery == 0
        @test i.force_stepfail
    end

    @testset "Filtered recovery preserves stage and estimator semantics" begin
        decay!(du, u, p, t) = (du .= -u; nothing)
        # Emulate a high-degree linearization that produces tiny corrections
        # without solving the stage equation. The value kernels remain exact;
        # the independent check must first try the usual endpoint refresh.
        for make in (TaylorGaussA4, ImplicitTaylor2A4), scale in (1.0e-9, 1.0, 1.0e6),
                linsolve in (nothing, ODTI.LinearSolve.SimpleLUFactorization())
            u0 = [scale]
            p = ODEProblem(decay!, u0, (0.0, 1.0))
            # The increment form: its Newton matrix comes from polynomial_jac.
            i = init(
                p, make(; linsolve, extrapolant = :constant, stages = :reduced);
                dt = 0.1, abstol = [0.01scale], reltol = 0.0
            )
            c = i.cache
            original = c.polynomial_jac
            # The replacement below deliberately violates the generated
            # derivative values, so its original constant-cache proof cannot
            # authorize reuse by the independent validation kernel.
            c.constant_jacobian = false
            c.polynomial_jac = typeof(original)(
                (b, b1, Jk, u, t, h) -> begin
                    original(b, b1, Jk, u, t, h)
                    Jk .*= 1.0e8
                    nothing
                end
            )
            @test ODTI.implicit_taylor_step!(i, c, 0.0, 0.1, u0, i.u)
            @test c.nendpoint_refresh == c.nfilter_recovery == 1
            @test c.nrestart == 0 && c.nvalidation_fail == 2
            @test c.ncoefjac + c.ncheckjac == i.stats.njacs
            @test c.iter == i.stats.nnonliniter <= i.alg.maxiters
            @test c.recovery_since_accept && c.newton_limit
            @test !c.jac_valid && !c.lu_valid
            @test c.filter_valid && !c.newton_failed && iszero(c.dt_cap)
            # Do not report a contraction rate measured before the direction
            # changed as a fast rate for the filtered correction.
            @test c.ηold == 1
            reference = if make === TaylorGaussA4
                scale * (1 - 0.1 / 2 + 0.1^2 / 12) / (1 + 0.1 / 2 + 0.1^2 / 12)
            else
                scale / (1 + 0.1c.μ + (0.1c.μ)^2 / 2)
            end
            @test abs(only(c.uintermediate) - reference) / (0.01scale) < c.κ
            @test c.embedded_filtered == (c.linsolve_filter.b isa AbstractMatrix)

            # Rebuilding the estimator inputs at the final stage must give
            # the same estimate as the cached companion / mu-Taylor term.
            saved = OrdinaryDiffEqCore.get_EEst(i)
            distance = make === TaylorGaussA4 ? 0.1 : (1 - c.μ) * 0.1
            tc = make === TaylorGaussA4 ? 0.1 : c.μ * 0.1
            c.polynomial_iter(c.tmp, c.tmp_embedded, c.uintermediate, tc, distance)
            c.embedded_filtered = false
            ODTI.implicit_taylor_error!(i, c, 0.0, 0.1, u0, i.u)
            @test OrdinaryDiffEqCore.get_EEst(i) ≈ saved rtol = 1.0e-12
        end
    end

    @testset "Filtered trials freeze the merit, consume budget, and roll back" begin
        stationary!(du, u, p, t) = (fill!(du, zero(eltype(du))); nothing)
        p = ODEProblem(stationary!, [0.0], (0.0, 1.0))
        weights = ODTI.pade_increment_coefficients(Val(2), Val(2))[2]
        # For f=0 the exact stage residual is z and its root is zero.
        # Give the kept filter an approximate inverse B. B=2 overshoots
        # with a full correction and reaches the root with half a correction.
        function trial_setup(B, z)
            i = init(p, TaylorGaussA4(); dt = 0.1, abstol = 1.0, reltol = 0.0)
            c = i.cache
            c.utilde .= z
            c.uintermediate .= z
            c.rhs .= 0
            c.rhs_embedded .= 0
            c.jac_valid = c.lu_valid = true
            c.stage_error .= z
            γ = weights[2] / 2 * 0.1
            filter_jac = fill((inv(sqrt(B)) - 1) / γ, 1, 1)
            ODTI.filter_error!(i, c, c.stage_error, 2, weights[2], 0.1, filter_jac)
            return i, abs(only(c.stage_error))
        end
        for budget in (0, 1, 2)
            i, defect = trial_setup(2.0, 1.0)
            c = i.cache
            nw, njacs, niter = i.stats.nw, i.stats.njacs, i.stats.nnonliniter
            ok, trials = ODTI.recover_filtered_stage!(
                i, c, 0.0, 0.1, 0.1, [0.0], weights, defect, budget
            )
            @test ok == (budget == 2)
            @test trials == budget
            @test i.stats.nnonliniter - niter == trials
            @test c.nfilter_recovery == (budget > 0 ? 1 : 0)
            if ok
                @test abs(only(c.utilde)) < c.κ
                @test c.uintermediate == c.utilde
                @test c.embedded_filtered
                @test i.stats.njacs == njacs # f=0: the valid constant derivative is reused
            else
                @test c.utilde == c.uintermediate == [1.0]
                # No Jacobian or factorization is needed while trying the
                # fixed direction. Exhaustion cannot trigger a fresh restart.
                @test i.stats.nw == nw && i.stats.njacs == njacs
                @test c.nrestart == 0
                if budget > 0
                    @test !c.jac_valid && !c.lu_valid && !c.filter_valid
                    @test !c.embedded_filtered
                end
            end
        end

        # B=0.6 produces a trial with frozen-filter defect 0.0096 < kappa,
        # but its true residual is 0.016. Refreshing the independent filter
        # must veto acceptance and restore the start, not trust the merit.
        i, defect = trial_setup(0.6, 0.04)
        c = i.cache
        ok, trials = ODTI.recover_filtered_stage!(
            i, c, 0.0, 0.1, 0.1, [0.0], weights, defect, 5
        )
        @test !ok && trials == 1
        @test c.ncheckjac == 0
        @test c.nvalidation_fail == 1
        @test c.utilde == c.uintermediate == [0.04]
        @test !c.jac_valid && !c.lu_valid && !c.filter_valid && !c.embedded_filtered
    end

    @testset "Filtered recovery cannot extend an exhausted attempt" begin
        decay!(du, u, p, t) = (du .= -u; nothing)
        p = ODEProblem(decay!, [1.0], (0.0, 1.0))
        for maxiters in (1, 2)
            i = init(
                p, TaylorGaussA4(; maxiters, extrapolant = :constant, stages = :reduced);
                dt = 0.1, abstol = 0.01, reltol = 0.0
            )
            c = i.cache
            original = c.polynomial_jac
            # As in the recovery test above, this fault-injected provider
            # no longer carries the generated constant-derivative proof.
            c.constant_jacobian = false
            c.polynomial_jac = typeof(original)(
                (b, b1, Jk, u, t, h) -> begin
                    original(b, b1, Jk, u, t, h)
                    Jk .*= 1.0e8
                    nothing
                end
            )
            @test !ODTI.implicit_taylor_stage!(i, c, 0.0, 0.1, [1.0], i.u)
            @test c.iter == i.stats.nnonliniter == maxiters
            @test c.nfilter_recovery == c.nrestart == 0
            @test i.force_stepfail && c.newton_failed
        end
    end

    @testset "Recovery feedback survives an error rejection" begin
        decay!(du, u, p, t) = (du .= -10 .* u; nothing)
        for controller in (:PI, :Predictive), direction in (-1, 1), inplace in (false, true)
            h = direction * 0.05
            tspan = (0.0, Float64(direction))
            p = inplace ? ODEProblem(decay!, [1.0], tspan) :
                ODEProblem((u, p, t) -> -10u, 1.0, tspan)
            # Feedback of the increment form's recovery paths.
            i = init(p, TaylorRadauLA3(; controller, stages = :reduced); dt = h, abstol = 1.0e-8, reltol = 1.0e-6)
            c = ODTI.newton_state(i.cache)
            # The same flag is set by both recovery routes. A subsequent
            # error-rejected attempt must not erase it when theta is reset.
            c.recovery_since_accept = true
            @test ODTI.implicit_taylor_step!(i, c, 0.0, h, [1.0], c.u)
            @test OrdinaryDiffEqCore.get_EEst(i) > 1
            @test c.recovery_since_accept && c.newton_limit
            @test ODTI.implicit_taylor_stage!(i, c, 0.0, h / 10, [1.0], c.u)
            @test c.recovery_since_accept && c.newton_limit
            c.dt_cap = abs(h)
            c.cap_probe, c.cap_streak = 2, 1
            OrdinaryDiffEqCore.set_EEst!(i, 0.01)
            next_h = OrdinaryDiffEqCore.step_accept_controller!(i, i.alg, 0.2)
            @test sign(next_h) == direction
            @test abs(next_h) <= abs(h)
            @test !c.recovery_since_accept && !c.newton_failed
            @test c.dt_cap == abs(h) && c.cap_streak == 0 && !c.cap_raised
            # A genuinely fast next stage may grow again.
            @test ODTI.implicit_taylor_stage!(i, c, 0.0, h / 10, [1.0], c.u)
            @test !c.newton_limit
            c.dt_cap = 0.0
            @test ODTI.newton_feedback!(c, 0.2, h) == 0.2

            # Reinitialization clears histories, diagnostics and feedback,
            # including the stage cache behind an out-of-place wrapper.
            c.recovery_since_accept = c.newton_failed = c.newton_limit = true
            c.newton_rate_measured = c.newton_slow = c.stage_handover = true
            c.dt_cap = 0.1
            c.nrestart = 3
            c.nfilter_recovery = 2
            SciMLBase.reinit!(i; reset_dt = false)
            @test !c.recovery_since_accept && !c.newton_failed && !c.newton_limit
            @test !c.newton_rate_measured && !c.newton_slow && !c.stage_handover
            @test !c.jac_valid && !c.lu_valid && !c.filter_valid
            @test iszero(c.dt_cap) && c.cap_probe == 1 && c.cap_streak == 0
            @test c.ncoefjac == c.ncheckjac == c.nrestart == c.nendpoint_refresh == c.nvalidation_fail == c.nfilter_recovery == 0
            step!(i)
            fresh = init(p, TaylorRadauLA3(; controller, stages = :reduced); dt = h, abstol = 1.0e-8, reltol = 1.0e-6)
            step!(fresh)
            @test i.t ≈ fresh.t
            @test i.u ≈ fresh.u
        end
    end
end
