using LinearAlgebra, ForwardDiff, OrdinaryDiffEqTaylorSeries

@testset "Implicit Taylor derivative stages" begin
    ODTI = OrdinaryDiffEqTaylorSeriesImplicit
    ODTS = OrdinaryDiffEqTaylorSeries
    k1, k2, k3 = 0.04, 3.0e7, 1.0e4
    rober!(du, u, p, t) = (
        du[1] = -k1 * u[1] + k3 * u[2] * u[3];
        du[2] = k1 * u[1] - k3 * u[2] * u[3] - k2 * u[2]^2; du[3] = k2 * u[2]^2; nothing
    )
    frober(u) = [-k1 * u[1] + k3 * u[2] * u[3], k1 * u[1] - k3 * u[2] * u[3] - k2 * u[2]^2, k2 * u[2]^2]
    # Taylor coefficients of f along the jet u + v1 σ + v2 σ²: F_0 = f(u), F_1 = f'(u) v1,
    # F_2 = f'(u) v2 + ½ f''(u)[v1, v1]; and J_m = [f_u]_m.
    function jet_coefficients(u, v1, v2)
        J = ForwardDiff.jacobian(frober, u)
        H(v) = ForwardDiff.jacobian(w -> ForwardDiff.jacobian(frober, w) * v, u)   # f''(u)[v, ·]
        F0 = frober(u); F1 = J * v1; F2 = J * v2 + H(v1) * v1 / 2
        return hcat(F0, F1, F2), cat(J, H(v1), H(v2); dims = 3)
    end

    @testset "Stage kernels: f and f_u along a free jet" begin
        values, jacobians = ODTS.build_stage_bundles(rober!, Val(true), nothing, 3, 3)
        rng = Random.Xoshiro(7)
        for _ in 1:5
            u = abs.(randn(rng, 3)); V = randn(rng, 3, 3) .* [1.0 1.0e-4 1.0e-8]
            F = zeros(3, 3); J = zeros(3, 3, 3); F2 = zeros(3, 3)
            jacobians(F, J, u, V, 0.0, 0.5)
            values(F2, u, V, 0.0, 0.5)
            Fe, Je = jet_coefficients(u, V[:, 1], V[:, 2])
            @test F ≈ Fe rtol = 1.0e-12
            @test F2 ≈ Fe rtol = 1.0e-12
            @test J ≈ Je rtol = 1.0e-12
        end
        @test !jacobians.constant_jacobian
        # A linear f has constant jet Jacobians: J_0 = A, J_1 = J_2 = 0.
        A = [-1.0 0.3; 0.0 -4.0]
        linear!(du, u, p, t) = (mul!(du, A, u); nothing)
        _, jl = ODTS.build_stage_bundles(linear!, Val(true), nothing, 3, 2)
        F = zeros(2, 3); J = zeros(2, 2, 3)
        jl(F, J, [1.0, 2.0], randn(rng, 2, 3), 0.0, 0.1)
        @test jl.constant_jacobian
        @test J[:, :, 1] ≈ A && iszero(J[:, :, 2]) && iszero(J[:, :, 3])
    end

    @testset "Requested stage Jacobian order preserves every jet value" begin
        calls = Ref(0)
        function nonlinear!(du, u, p, t)
            calls[] += 1
            du[1] = p[1] * u[1]^2 + t * u[2]
            du[2] = u[1] * u[2] + p[2] * t^2
            return nothing
        end
        function nonlinear(u, p, t)
            calls[] += 1
            return [p[1] * u[1]^2 + t * u[2], u[1] * u[2] + p[2] * t^2]
        end
        # Direct polynomial convolution, independent of generated differentiation.
        function expected_jet(u, V, p, t, h)
            Q = size(V, 2)
            U = hcat(u, V[:, 1:(Q - 1)])
            F = zeros(2, Q)
            J = zeros(2, 2, Q)
            for k in 0:(Q - 1)
                F[1, k + 1] = p[1] * sum(U[1, j + 1] * U[1, k - j + 1] for j in 0:k) + t * U[2, k + 1]
                k >= 1 && (F[1, k + 1] += h * U[2, k])
                F[2, k + 1] = sum(U[1, j + 1] * U[2, k - j + 1] for j in 0:k) +
                    p[2] * (k == 0 ? t^2 : k == 1 ? 2t * h : k == 2 ? h^2 : 0.0)
                J[1, 1, k + 1] = 2p[1] * U[1, k + 1]
                J[1, 2, k + 1] = k == 0 ? t : k == 1 ? h : 0.0
                J[2, 1, k + 1] = U[2, k + 1]
                J[2, 2, k + 1] = U[1, k + 1]
            end
            return F, J
        end
        p, other_p = (1.7, -0.4), (-0.6, 2.1)
        u, t, h = [0.8, -0.3], 0.27, -0.15
        for iip in (true, false), Q in (1, 2, 3, 4)
            f = iip ? nonlinear! : nonlinear
            V = reshape(collect(range(-0.17, 0.23; length = 2Q)), 2, Q)
            full = ODTS.build_stage_bundles(f, Val(iip), p, Q, 2)
            only0 = ODTS.build_stage_bundles(f, Val(iip), p, Q, 2; jacobian_order = 0)
            @test full === ODTS.build_stage_bundles(f, Val(iip), p, Q, 2; jacobian_order = Q - 1)
            @test (full === only0) == (Q == 1)
            # Prewarm the ODEFunction entry separately: unwrapped_f preserves its
            # wrapper, so its cache key need not be the same as the bare callable.
            wrappedf = SciMLBase.ODEFunction{iip}(f)
            wrapped0 = ODTS.build_stage_bundles(wrappedf, p, Q, 2; jacobian_order = 0)
            # Repeating either already-built entry must not regenerate symbols.
            previous_calls = calls[]
            @test only0 === ODTS.build_stage_bundles(f, Val(iip), p, Q, 2; jacobian_order = 0)
            @test wrapped0 === ODTS.build_stage_bundles(wrappedf, p, Q, 2; jacobian_order = 0)
            @test calls[] == previous_calls
            F, J = zeros(2, Q), zeros(2, 2, Q)
            F0, J0 = zeros(2, Q), zeros(2, 2, 1)
            Fvalues = zeros(2, Q)
            full[2](F, J, u, V, t, h)
            only0[2](F0, J0, u, V, t, h)
            only0[1](Fvalues, u, V, t, h)
            expected_F, expected_J = expected_jet(u, V, p, t, h)
            @test F ≈ expected_F rtol = 1.0e-12 atol = 1.0e-14
            @test F0 ≈ expected_F rtol = 1.0e-12 atol = 1.0e-14
            @test Fvalues ≈ expected_F rtol = 1.0e-12 atol = 1.0e-14
            @test J ≈ expected_J rtol = 1.0e-12 atol = 1.0e-14
            @test J0[:, :, 1] ≈ expected_J[:, :, 1] rtol = 1.0e-12
            @test !full[2].constant_jacobian && !only0[2].constant_jacobian
            # The last free stage is absent from every coefficient of f used here.
            V[:, end] .+= 10
            only0[2](F0, J0, u, V, t, h)
            @test F0 ≈ expected_F rtol = 1.0e-12 atol = 1.0e-14
            @test J0[:, :, 1] ≈ expected_J[:, :, 1] rtol = 1.0e-12
            if Q == 3
                alternate = ODTS.build_stage_bundles(f, Val(iip), other_p, Q, 2; jacobian_order = 0)
                @test alternate !== only0
                alternate[2](F0, J0, u, V, t, h)
                alternate_F, alternate_J = expected_jet(u, V, other_p, t, h)
                @test F0 ≈ alternate_F rtol = 1.0e-12 atol = 1.0e-14
                @test J0[:, :, 1] ≈ alternate_J[:, :, 1] rtol = 1.0e-12
                only0[2](F0, J0, u, V, t, h)
                @test F0 ≈ expected_F rtol = 1.0e-12 atol = 1.0e-14
                @test J0[:, :, 1] ≈ expected_J[:, :, 1] rtol = 1.0e-12
            end
            if !iip && Q == 4
                _, first_two = ODTS.build_stage_bundles(f, Val(false), p, Q, 2; jacobian_order = 1)
                J01 = zeros(2, 2, 2)
                first_two(F0, J01, u, V, t, h)
                @test F0 ≈ expected_F rtol = 1.0e-12 atol = 1.0e-14
                @test J01 ≈ expected_J[:, :, 1:2] rtol = 1.0e-12 atol = 1.0e-14
            end
            previous_calls = calls[]
            @test_throws ArgumentError ODTS.build_stage_bundles(f, Val(iip), p, Q, 2; jacobian_order = -1)
            @test_throws ArgumentError ODTS.build_stage_bundles(f, Val(iip), p, Q, 2; jacobian_order = Q)
            @test calls[] == previous_calls
        end
        @test_throws ArgumentError ODTS.build_stage_bundles(nonlinear!, Val(true), p, 0, 2)
        # Time-dependent forcing does not invalidate a constant state Jacobian.
        affine(u, p, t) = [p[1] * u[1] + u[2] + t^2, -u[2] + t]
        varying(u, p, t) = [t * u[1] + u[2], u[1] - t * u[2]]
        for Q in (1, 3), order in unique((0, Q - 1))
            _, constant = ODTS.build_stage_bundles(affine, Val(false), p, Q, 2; jacobian_order = order)
            _, nonconstant = ODTS.build_stage_bundles(varying, Val(false), nothing, Q, 2; jacobian_order = order)
            @test constant.constant_jacobian
            @test !nonconstant.constant_jacobian
            F, J = zeros(2, Q), zeros(2, 2, order + 1)
            constant(F, J, u, zeros(2, Q), t, h)
            @test J[:, :, 1] ≈ [p[1] 1.0; 0.0 -1.0]
            @test all(iszero, J[:, :, 2:end])
        end
    end

    @testset "Simplified Newton matrix splits on the Padé roots" begin
        rng = Random.Xoshiro(8)
        A = [-1.0 0.3 0.0; 0.0 -4.0 1.0; 0.5 0.0 -20.0]
        linear!(du, u, p, t) = (mul!(du, A, u); nothing)
        integrator = init(ODEProblem(linear!, ones(3), (0.0, 1.0)), TaylorRadauLA5())
        @test size(integrator.cache.Jst) == (3, 3, 1)
        for (P, Q) in ((2, 2), (1, 3), (2, 3), (2, 4), (3, 5))
            b = ODTI.pade_increment_coefficients(Val(P), Val(Q))[2]
            C = zeros(Q, Q); C[1, :] .= b[2:end]
            for k in 2:Q
                C[k, k - 1] = -1 / k
            end
            # The eigenvalues of C are -1/r over the roots r of the Padé denominator.
            λr, λc = ODTI.stage_eigensystem(b)
            roots = ODTI.denominator_roots(P, Q)
            # Match each eigenvalue to its nearest -1/r (a sort can split a conjugate pair).
            λ = vcat(λr, λc, conj.(λc)); μ = -1 ./ roots
            @test length(λ) == length(μ) == Q
            @test all(x -> minimum(abs.(x .- μ)) <= 1.0e-10 * abs(x), λ)
            @test all(y -> minimum(abs.(y .- λ)) <= 1.0e-10 * abs(y), μ)
            # The split solve inverts I + C ⊗ hJ.
            n = 5; h = 0.37
            J = randn(rng, n, n) - 2I
            S = ODTI.build_stage_factors(b, J, nothing, zeros(n), integrator.p, integrator.t)
            ODTI.form_stage_factors!(S, J, h)
            r = randn(rng, Q * n); dv = zeros(Q * n)
            ODTI.solve_stage_factors!(integrator, S, dv, r, n, true)
            @test dv ≈ (I + kron(C, h * J)) \ r rtol = 1.0e-9
            # Kept factors: a second right-hand side without refactorizing.
            r2 = randn(rng, Q * n)
            ODTI.solve_stage_factors!(integrator, S, dv, r2, n, false)
            @test dv ≈ (I + kron(C, h * J)) \ r2 rtol = 1.0e-9
        end
    end

    @testset "Factors are kept while I + λ_i hJ changes little" begin
        # stage_factors_outdated against the row rule applied to every factor directly.
        rng = Random.Xoshiro(10)
        A = [-1.0 0.3 0.0; 0.0 -4.0 1.0; 0.5 0.0 -20.0]
        linear!(du, u, p, t) = (mul!(du, A, u); nothing)
        integrator = init(ODEProblem(linear!, ones(3), (0.0, 1.0)), TaylorRadauLA5())
        cutoff = 1 / 20
        for (P, Q) in ((2, 2), (1, 3), (2, 3))
            b = ODTI.pade_increment_coefficients(Val(P), Val(Q))[2]
            n = 6; h = 0.37
            J = randn(rng, n, n) - 5I
            S = ODTI.build_stage_factors(b, J, nothing, zeros(n), integrator.p, integrator.t)
            @test ODTI.stage_factors_outdated(S, J, h, cutoff)      # nothing factorized yet
            ODTI.form_stage_factors!(S, J, h)
            @test !ODTI.stage_factors_outdated(S, J, h, cutoff)
            λs = vcat(S.λr, S.λc, conj.(S.λc))
            outcomes = Bool[]
            for scale in (1.0e-4, 1.0e-3, 1.0e-2, 3.0e-2, 0.1, 0.3), _ in 1:4
                J2 = J .+ scale .* randn(rng, n, n)
                h2 = h * (1 + scale * randn(rng))
                expected = any(λs) do λ
                    E = λ .* (h2 .* J2 .- h .* J)
                    M = I + λ .* h2 .* J2
                    any(k -> sum(abs, E[k, :]) > cutoff * max(1, sum(abs, M[k, :])), 1:n)
                end
                @test ODTI.stage_factors_outdated(S, J2, h2, cutoff) == expected
                push!(outcomes, expected)
            end
            @test any(outcomes) && !all(outcomes)
        end
    end

    @testset "Exact Newton matrix for linear f" begin
        # f_u is constant, so the left-out blocks J_m (m ≥ 1) vanish: with z eliminated the
        # Jacobian of the block residual is I + C ⊗ hJ and the split solve is Newton's.
        rng = Random.Xoshiro(9)
        n, P, Q, h = 3, 2, 3, 0.7
        L = randn(rng, n, n)
        b = ODTI.pade_increment_coefficients(Val(P), Val(Q))[2]
        uprev = randn(rng, n); A = randn(rng, n)
        function residual(x)
            z = x[1:n]; V = reshape(x[(n + 1):end], n, Q)
            jet = [L * (uprev .+ z), (L * V[:, k] for k in 1:(Q - 1))...]
            R0 = z .+ sum(b[k + 1] .* V[:, k] for k in 1:Q) .- A
            vcat(R0, (V[:, k] .- (h / k) .* jet[k] for k in 1:Q)...)
        end
        x = randn(rng, (Q + 1) * n)
        JB = ForwardDiff.jacobian(residual, x)
        S = JB[(n + 1):end, (n + 1):end] - JB[(n + 1):end, 1:n] * JB[1:n, (n + 1):end]
        C = zeros(Q, Q); C[1, :] .= b[2:end]
        for k in 2:Q
            C[k, k - 1] = -1 / k
        end
        @test S ≈ I + kron(C, h * L) rtol = 1.0e-12
        # The right-hand side and the recovered δz give the full Newton correction.
        z = x[1:n]; V = reshape(x[(n + 1):end], n, Q)
        Fst = hcat(L * (uprev .+ z), (L * V[:, k] for k in 1:(Q - 1))...)
        rb = zeros(Q * n); r0 = zeros(n); dz = zeros(n)
        ODTI.stage_residual!(rb, r0, z, V, Fst, L, A, b, h)
        dv = S \ rb
        ODTI.stage_increment_correction!(dz, r0, dv, V, b)
        @test vcat(dz, dv) ≈ JB \ residual(x) rtol = 1.0e-9
    end

    @testset "Same Padé step as the increment form" begin
        A = [-1.0 0.3 0.0; 0.0 -4.0 1.0; 0.5 0.0 -20.0]
        linear!(du, u, p, t) = (mul!(du, A, u); nothing)
        u0 = [1.0, 0.5, -0.25]
        prob = ODEProblem(linear!, u0, (0.0, 0.5))
        for make in (TaylorRadauLA3, TaylorRadauLA5, TaylorLobattoLA4, TaylorLobattoLA6, TaylorGaussA4)
            kw = (dt = 0.05, adaptive = false)
            d = solve(prob, make(stages = :derivative); kw...)
            r = solve(prob, make(stages = :reduced); kw...)
            @test SciMLBase.successful_retcode(d)
            @test d.u[end] ≈ r.u[end] rtol = 1.0e-9
            @test norm(d.u[end] - exp(0.5A) * u0) < 1.0e-4
        end
        # Constant jet Jacobians are not re-evaluated per step.
        counts = map((0.05, 0.01)) do dt
            i = init(prob, TaylorRadauLA5(stages = :derivative); dt, adaptive = false)
            @test i.cache.stage_system
            @test i.cache.constant_stage_jacobian
            solve!(i)
            i.cache.ncoefjac
        end
        @test counts[1] == counts[2]
    end

    @testset "Accepted stages solve the increment form" begin
        # Compare both stage equations at identical starting states and step
        # sizes. dtmax keeps the steps within the reduced solve's convergence range.
        prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e5))
        i = init(
            prob, TaylorRadauLA5(stages = :derivative); abstol = 1.0e-6, reltol = 1.0e-3,
            dtmax = 0.25, save_everystep = false
        )
        n = 0
        while i.t < 1.0
            uprev = copy(i.u); t = i.t
            step!(i)
            h = i.cache.stages_h
            r = solve(
                remake(prob; u0 = uprev, tspan = (t, t + h)), TaylorRadauLA5(stages = :reduced);
                dt = h, adaptive = false, save_everystep = false
            )
            tol = 1.0e-6 .+ 1.0e-3 .* abs.(i.u)
            @test maximum(abs.(r.u[end] .- i.u) ./ tol) < 0.05
            n += 1
        end
        @test n >= 10
    end

    @testset "ROBER is limited by the error estimate, not by Newton" begin
        prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e5))
        for make in (TaylorRadauLA5, TaylorLobattoLA4)
            d = solve(prob, make(stages = :derivative); abstol = 1.0e-6, reltol = 1.0e-3, save_everystep = false)
            r = solve(prob, make(stages = :reduced); abstol = 1.0e-6, reltol = 1.0e-3, save_everystep = false)
            @test SciMLBase.successful_retcode(d)
            @test d.stats.naccept < 200
            @test 3 * d.stats.naccept < r.stats.naccept
            @test d.stats.nnonlinconvfail < 10
            @test d.u[end] ≈ r.u[end] atol = 2.0e-3
        end
    end

    @testset "Auto mode uses the accepted controller's growth proposal" begin
        decay!(du, u, p, t) = (du .= -u; nothing)
        # History can turn a raw growth request into shrinkage or a steady
        # step. The steady-step interval can suppress raw PI growth too.
        cases = (
            (:Predictive, 0.5, 1.0, 1.0, false),
            (:Predictive, 0.5, sqrt(0.01 * (0.9 * 1.1)^5), 1.0, false),
            (:Predictive, 0.5, 0.01, 1.0, true),
            (:Predictive, 0.95, 0.01, 0.9, false),
            (:PI, 0.95, 0.01, 0.9, false),
            (:PI, 0.5, 0.01, 1.0, true),
            (:PI, 1.5, 0.01, 1.0, false),
        )
        for (controller, q, eest, steady_min, grows) in cases,
                direction in (-1, 1), inplace in (false, true),
                cause in (:fast, :slow, :recovery, :recovered_slow)
            h = direction * 0.1
            p = inplace ? ODEProblem(decay!, [1.0], (0.0, Float64(direction))) :
                ODEProblem((u, p, t) -> -u, 1.0, (0.0, Float64(direction)))
            alg = TaylorRadauLA5(; controller)
            make() = begin
                ctl = controller === :PI ?
                    OrdinaryDiffEqCore.PIController(alg; qsteady_min = steady_min) :
                    OrdinaryDiffEqCore.PredictiveController(; qsteady_min = steady_min)
                i = init(p, alg; dt = h, controller = ctl)
                i.success_iter = 1
                OrdinaryDiffEqCore.set_EEst!(i, eest)
                if controller === :Predictive
                    i.controller_cache.dtacc = h
                    i.controller_cache.erracc = 0.01
                end
                i
            end
            reference = make()
            # The underlying controller is the oracle: it resolves its own
            # history and deadband, without Taylor's Newton growth limit.
            proposed = invoke(
                OrdinaryDiffEqCore.step_accept_controller!,
                Tuple{Any, typeof(reference.controller_cache), Any, Any},
                reference, reference.controller_cache, alg, q
            )
            @test (abs(proposed) > abs(h)) == grows

            i = make()
            c = ODTI.newton_state(i.cache)
            slow = cause in (:slow, :recovered_slow)
            c.recovery_since_accept = cause in (:recovery, :recovered_slow)
            ODTI.record_newton_convergence!(
                c, cause !== :recovery,
                slow || cause === :recovery ? 1.0 : 0.0, alg
            )
            next_h = OrdinaryDiffEqCore.step_accept_controller!(i, alg, q)
            held = cause !== :fast
            expected = held && grows ? h : proposed
            @test next_h == expected
            @test c.use_stages == (held && grows)
            @test c.newton_rate_measured == (cause !== :recovery)
            @test c.newton_slow == slow
            @test !c.recovery_since_accept && !c.newton_failed
            # Controller history must agree with one underlying invocation.
            if controller === :Predictive
                @test i.controller_cache.dtacc == reference.controller_cache.dtacc
                @test i.controller_cache.erracc == reference.controller_cache.erracc
            else
                @test i.controller_cache.errold == reference.controller_cache.errold
            end
        end
    end

    @testset "stages = :auto follows the increment form's Newton iteration" begin
        # Where the increment form converges, :auto is the increment form.
        A = [-1.0 0.3 0.0; 0.0 -4.0 1.0; 0.5 0.0 -20.0]
        linear!(du, u, p, t) = (mul!(du, A, u); nothing)
        prob = ODEProblem(linear!, [1.0, 0.5, -0.25], (0.0, 2.0))
        i = init(prob, TaylorRadauLA5(); abstol = 1.0e-8, reltol = 1.0e-6)
        c = ODTI.newton_state(i.cache)
        @test c.stage_system && c.stage_auto && !c.use_stages
        solve!(i)
        r = solve(prob, TaylorRadauLA5(stages = :reduced); abstol = 1.0e-8, reltol = 1.0e-6)
        @test !c.use_stages && c.stage_probe == 1
        @test i.sol.t == r.t && i.sol.u == r.u

        # A failed increment solve hands its step to the derivative stages at the same dt:
        # no dt cap, no failed step. From ROBER's state at t = 10 a step of 1 is beyond the
        # increment form's Newton basin but not the derivative stages'.
        rober = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e5))
        u10 = solve(
            rober, TaylorRadauLA5(stages = :derivative); abstol = 1.0e-12, reltol = 1.0e-10,
            save_everystep = false, saveat = [10.0]
        ).u[1]
        step(st; kw...) = begin
            j = init(
                remake(rober; u0 = u10, tspan = (10.0, 1.0e5)), TaylorRadauLA5(stages = st);
                dt = 1.0, abstol = 1.0e-8, reltol = 1.0e-5, kw...
            )
            out = similar(j.u)
            ok = ODTI.implicit_taylor_step!(j, j.cache, 10.0, 1.0, copy(j.u), out)
            j, ok, out
        end
        _, ok_r, _ = step(:reduced)
        _, ok_d, out_d = step(:derivative)
        @test !ok_r && ok_d
        j, ok, out = step(:auto)
        c = ODTI.newton_state(j.cache)
        @test ok && c.use_stages
        @test iszero(c.dt_cap) && !c.newton_failed
        @test j.stats.nnonlinconvfail == 1
        @test out ≈ out_d rtol = 1.0e-10
        # Fixed steps hand over that step only.
        j, ok, out = step(:auto; adaptive = false)
        c = ODTI.newton_state(j.cache)
        @test ok && !c.use_stages
        @test out ≈ step(:derivative; adaptive = false)[3] rtol = 1.0e-10

        # Mode bookkeeping (stage_mode_feedback!, q < 1 asks for growth).
        c = ODTI.newton_state(init(prob, TaylorRadauLA5()).cache)
        fast!() = begin
            c.newton_limit = c.newton_slow = c.newton_failed = c.recovery_since_accept = false
            c.newton_rate_measured = true
        end
        ODTI.hand_over_to_stages!(c, true)
        @test c.use_stages && c.stage_probe == 1          # not a trial: no backoff
        fast!(); ODTI.stage_mode_feedback!(c, 0.5)
        @test !c.use_stages && c.stage_trial && !c.stage_handover # accepted rescue counts toward the probe
        ODTI.hand_over_to_stages!(c, true)                 # the trial fails
        @test c.use_stages && !c.stage_trial && c.stage_probe == 2
        fast!(); ODTI.stage_mode_feedback!(c, 0.5)
        @test c.use_stages && c.stage_streak == 1 && !c.stage_handover # rescue counts once
        c.newton_limit = c.newton_slow = true; ODTI.stage_mode_feedback!(c, 0.5)   # a slow step restarts the count
        @test c.use_stages && c.stage_streak == 0
        fast!(); ODTI.stage_mode_feedback!(c, 0.5); ODTI.stage_mode_feedback!(c, 0.5)
        @test !c.use_stages && c.stage_trial
        c.newton_limit = c.newton_slow = true; ODTI.stage_mode_feedback!(c, 0.5)   # the trial holds back growth
        @test c.use_stages && c.stage_probe == 4
        fast!()
        for _ in 1:4
            ODTI.stage_mode_feedback!(c, 0.5)
        end
        @test !c.use_stages && c.stage_trial
        ODTI.stage_mode_feedback!(c, 0.5)                  # the trial converges fast
        @test !c.use_stages && !c.stage_trial && c.stage_probe == 1
        c.newton_limit = c.newton_slow = true; ODTI.stage_mode_feedback!(c, 1.5)   # slow, but no growth wanted
        @test !c.use_stages
        ODTI.stage_mode_feedback!(c, 0.5)                  # slow and holding back growth
        @test c.use_stages && c.stage_probe == 1

        # A rescue survives error retries until acceptance, but is then consumed.
        ODTI.hand_over_to_stages!(c, true)
        @test c.stage_handover
        c.recovery_since_accept = true
        fast_rate = 0.0
        ODTI.record_newton_convergence!(c, true, fast_rate, TaylorRadauLA5())
        @test c.stage_handover && c.newton_limit && !c.newton_slow
        ODTI.record_newton_convergence!(c, false, 1.0, TaylorRadauLA5())
        @test c.stage_handover && !c.newton_rate_measured && !c.newton_slow
        ODTI.newton_feedback!(c, 0.5, 0.1)
        @test c.use_stages && !c.stage_trial && !c.stage_handover && c.stage_streak == 0
        @test !c.recovery_since_accept
        fast!(); ODTI.stage_mode_feedback!(c, 0.5)
        @test !c.use_stages && c.stage_trial

        # Recovery remains a switching signal even if the final rate is fast or unmeasured.
        for measured in (false, true)
            c.use_stages = c.stage_trial = false
            c.recovery_since_accept = true
            ODTI.record_newton_convergence!(c, measured, measured ? 0.0 : 1.0, TaylorRadauLA5())
            @test c.newton_limit && !c.newton_slow
            @test ODTI.newton_feedback!(c, 0.5, 0.1) == 1.0
            @test c.use_stages && !c.newton_slow
        end

        # ROBER: the derivative stages take over, trials back off.
        a = solve(rober, TaylorRadauLA5(); abstol = 1.0e-6, reltol = 1.0e-3, save_everystep = false)
        r = solve(rober, TaylorRadauLA5(stages = :reduced); abstol = 1.0e-6, reltol = 1.0e-3, save_everystep = false)
        @test SciMLBase.successful_retcode(a)
        @test 3 * a.stats.naccept < r.stats.naccept
        k = init(rober, TaylorRadauLA5(); abstol = 1.0e-6, reltol = 1.0e-3, save_everystep = false)
        solve!(k)
        c = ODTI.newton_state(k.cache)
        @test c.use_stages && c.stage_probe >= 8
        SciMLBase.reinit!(k)                               # reinit starts over on the increment form
        @test !c.use_stages && !c.stage_trial && c.stage_probe == 1 && c.stage_streak == 0
        @test !c.stage_handover && !c.newton_rate_measured && !c.newton_slow
    end

    @testset "Stiff Newton error leaks into the next step's explicit side" begin
        # The estimate of stage_leak_norm! against the same formula on dense matrices:
        # W₁⁻¹ Σ_{k≤P} a_k h^k/k! [(J+ΔJ)^k - J^k] δ, W₁ = I + (b₁/Q) h J, in the tolerance norm.
        rng = Random.Xoshiro(11)
        prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e5))
        for make in (TaylorRadauLA5, TaylorLobattoLA4, TaylorGaussA4)
            alg = make(stages = :derivative)
            P, Q = ODTI._unwrap_val(alg.order), ODTI._unwrap_val(alg.order_q)
            a, b = ODTI.pade_increment_coefficients(Val(P), Val(Q))
            i = init(prob, alg; dt = 0.1, abstol = 1.0e-8, reltol = 1.0e-5)
            c = i.cache
            n = 3; h = 37.0
            J = randn(rng, n, n) - 50I; Jprev = J .+ 0.1 .* randn(rng, n, n)
            δ = 1.0e-9 .* randn(rng, n)
            Jst = reshape(copy(J), n, n, 1); Jp = reshape(copy(Jprev), n, n, 1)
            uprev = [0.5, 1.0e-6, 0.5]; u = uprev .+ 1.0e-3
            leak = ODTI.stage_leak_norm!(
                i, c, δ, @view(Jst[:, :, 1]), @view(Jp[:, :, 1]), h, Q, a, b[2],
                uprev, u, i.opts.abstol, i.opts.reltol, i.opts.internalnorm, 0.0
            )
            ΔJ = J - Jprev
            ΔE = sum(a[k + 1] * h^k / factorial(k) * ((J + ΔJ)^k - J^k) * δ for k in 1:P)
            expected = (I + (b[2] / Q) * h * J) \ ΔE
            scale = i.opts.abstol .+ i.opts.reltol .* max.(abs.(uprev), abs.(u))
            @test leak ≈ i.opts.internalnorm(expected ./ scale, 0.0) rtol = 1.0e-10
            @test c.filter_valid                     # the estimate's filter was factorized here
        end

        # Compare the predicted leak with the numerical effect of perturbing
        # ROBER's stiff component, using the reference flow as a control.
        ref = solve(prob, Rodas5P(); abstol = 1.0e-15, reltol = 1.0e-13)
        rt = 1.0e-4; at = rt / 1000
        t0 = 6218.0; h = 1459.0
        u0 = ref(t0)
        exact(x) = solve(
            remake(prob; u0 = x, tspan = (t0, t0 + h)), Rodas5P(); abstol = 1.0e-16, reltol = 1.0e-14,
            save_everystep = false
        ).u[end]
        function pade_step(x)
            j = init(
                remake(prob; u0 = x, tspan = (t0, 1.0e5)), TaylorRadauLA5(stages = :derivative, κ = 1.0e-11, maxiters = 300);
                dt = h, abstol = at, reltol = rt
            )
            out = similar(j.u)
            @test ODTI.implicit_taylor_step!(j, j.cache, t0, h, copy(x), out)
            j, out
        end
        _, p0 = pade_step(u0)
        δ = 3.0e-14
        u1 = copy(u0); u1[2] += δ
        j1, p1 = pade_step(u1)
        scale = at .+ rt .* abs.(u0)
        @test maximum(abs.(exact(u1) .- exact(u0)) ./ scale) < 1.0e-3
        measured = maximum(abs.((p1 .- p0)) ./ scale)
        @test measured > 1.0                                   # several tol from 3e-7 tol
        c = j1.cache
        a, b = ODTI.pade_increment_coefficients(Val(2), Val(3))
        F = similar(c.Fst); Jn = similar(c.Jst); Jp = similar(c.Jst)
        c.stage_jac(F, Jn, p1, c.V, t0 + h, h)
        c.stage_jac(F, Jp, u1, c.V, t0, h)
        e = zeros(3); e[2] = δ
        estimate = ODTI.stage_leak_norm!(
            j1, c, e, @view(Jn[:, :, 1]), @view(Jp[:, :, 1]), h, 3, a, b[2],
            u1, p1, at, rt, j1.opts.internalnorm, t0
        )
        @test 0.1 < estimate / measured < 10

        # Check the endpoint error and rejection count with leak control enabled.
        s = solve(prob, TaylorRadauLA5(); abstol = at, reltol = rt, save_everystep = false)
        err = maximum(abs.(s.u[end] .- ref.u[end]) ./ (at .+ rt .* abs.(ref.u[end])))
        @test err < 0.5
        @test ODTI.newton_state(init(prob, TaylorRadauLA5(); abstol = at, reltol = rt).cache).nleak == 0
        k = init(prob, TaylorRadauLA5(); abstol = at, reltol = rt, save_everystep = false)
        solve!(k)
        @test ODTI.newton_state(k.cache).nleak > 10                # the test does bind on ROBER
        s8 = solve(prob, TaylorRadauLA5(); abstol = 1.0e-11, reltol = 1.0e-8, save_everystep = false)
        @test s8.stats.nreject < 10
        # Constant Jacobians skip the test.
        A = [-1.0 0.3; 0.0 -4.0]
        linear!(du, u, p, t) = (mul!(du, A, u); nothing)
        l = init(ODEProblem(linear!, [1.0, 0.5], (0.0, 1.0)), TaylorRadauLA5(stages = :derivative); abstol = 1.0e-9, reltol = 1.0e-6)
        solve!(l)
        @test l.cache.constant_stage_jacobian && l.cache.nleak == 0 && !l.cache.jprev_valid
    end

    @testset "reinit and discontinuities drop the stages" begin
        prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 100.0))
        i = init(prob, TaylorRadauLA5(stages = :derivative); abstol = 1.0e-8, reltol = 1.0e-6)
        c = i.cache
        step!(i); step!(i)
        @test c.stages_valid && c.lu_stages_valid
        SciMLBase.derivative_discontinuity!(i, true)
        step!(i)
        @test c.stages_valid
        SciMLBase.reinit!(i)
        @test !c.stages_valid && !c.lu_stages_valid && !c.stage_jac_valid
        solve!(i)
        @test SciMLBase.successful_retcode(i.sol)
    end
end
