@testset "Implicit Taylor cross-step Jacobian reuse" begin
    ODTI = OrdinaryDiffEqTaylorSeriesImplicit

    @testset "Constant linearization is reused across accepted steps" begin
        A = [-1.0 0.3; 0.0 -4.0]
        linear!(du, u, p, t) = (mul!(du, A, u); nothing)
        u0 = [1.0, 0.5]
        for make in (
                TaylorGaussA4, TaylorRadauLA3, TaylorRadauLA5,
                TaylorLobattoLA4, ImplicitTaylor1A2, ImplicitTaylor2A4,
            )
            p = ODEProblem(linear!, u0, (0.0, 0.1))
            i = init(p, make(); dt = 0.01, dtmax = 0.01, abstol = 1.0e-9, reltol = 1.0e-6)
            solve!(i)
            c = ODTI.newton_state(i.cache)
            @test SciMLBase.successful_retcode(i.sol)
            @test i.stats.naccept >= 10
            if c.stage_system && c.use_stages
                # Derivative stages: the jet Jacobians are constant and evaluated once.
                @test c.constant_stage_jacobian
                @test c.ncoefjac == 1
            else
                @test c.constant_jacobian
                @test c.ncoefjac == 1
                @test c.ncheckjac == 0
            end
            @test i.u ≈ exp(0.1A) * u0 rtol = 2.0e-5
            @test !c.newton_failed && iszero(c.dt_cap)
        end
    end

    @testset "Variable derivatives retain fresh-step linearizations" begin
        square!(du, u, p, t) = (du .= u .^ 2; nothing)
        p = ODEProblem(square!, [0.1], (0.0, 1.0))
        i = init(p, TaylorGaussA4(); dt = 0.01, abstol = 1.0e-8, reltol = 1.0e-6)
        c = i.cache
        @test !c.constant_jacobian
        @test ODTI.implicit_taylor_stage!(i, c, 0.0, 0.01, [0.1], i.u)
        first = copy(i.u)
        ncoef = c.ncoefjac
        @test ODTI.implicit_taylor_stage!(i, c, 0.01, 0.01, first, i.u)
        @test c.ncoefjac == ncoef + 1
        @test c.ncheckjac >= 2
    end

    @testset "Time-dependent forcing can have constant coefficient Jacobians" begin
        affine!(du, u, p, t) = (du .= -u .+ t; nothing)
        p = ODEProblem(affine!, [1.0], (0.0, 1.0))
        i = init(p, TaylorRadauLA5(); dt = 0.02, dtmax = 0.02, abstol = 1.0e-10, reltol = 1.0e-8)
        solve!(i)
        @test i.cache.constant_jacobian
        @test i.cache.ncoefjac == 1
        @test i.cache.ncheckjac == 0
        @test maximum(abs(only(u) - (t - 1 + 2exp(-t))) for (u, t) in zip(i.sol.u, i.sol.t)) < 1.0e-8
    end

    @testset "Constant-derivative checks retain values, filtering, and invalidation" begin
        # The derivative is constant, but the forcing and both polynomials
        # change with the evaluation time and distance. Compare with an
        # independently evaluated endpoint Jacobian, including degree one
        # and complex mu-Taylor expansion centers.
        affine!(du, u, p, t) = (@. du = p[1] * u + p[2] * sin(t); nothing)
        uprev = [0.5, -0.1]
        p = ODEProblem(affine!, uprev, (0.0, 1.0), [-2.0, 0.3])
        for make in (ImplicitTaylor1L1, TaylorGaussA4, TaylorRadauLA5, ImplicitTaylor2A4)
            alg = make(stages = :reduced)
            i = init(p, alg; dt = 0.02, abstol = 1.0e-9, reltol = 1.0e-6)
            reference = init(p, alg; dt = 0.02, abstol = 1.0e-9, reltol = 1.0e-6)
            for j in (i, reference)
                @test ODTI.implicit_taylor_step!(j, j.cache, 0.0, 0.02, uprev, j.u)
                @test j.cache.constant_jacobian && j.cache.jac_valid
            end
            reference.cache.constant_jacobian = false
            c, cr = i.cache, reference.cache
            stored = copy(c.Jk)
            P, Q = ODTI._unwrap_val(alg.order), ODTI._unwrap_val(alg.order_q)
            pade = ODTI.is_taylor_pade(alg)
            weights = pade ? ODTI.pade_increment_coefficients(Val(P), Val(Q))[2] :
                ODTI.taylor_increment_coefficients(Val(P))
            for h in (0.02, 0.035)
                tc = pade ? 0.1 + h : 0.1 + c.μ * h
                distance = pade ? h : -c.μ * h
                for cc in (c, cr)
                    cc.utilde .= [0.02, -0.03]
                    cc.uintermediate .= uprev .+ cc.utilde
                end
                njacs, ncheck = i.stats.njacs, c.ncheckjac
                expected = ODTI.check_stage_defect!(reference, cr, 0.1, tc, distance, uprev, weights)
                observed = ODTI.check_stage_defect!(i, c, 0.1, tc, distance, uprev, weights)
                @test observed[1] ≈ expected[1] rtol = 1.0e-12
                @test observed[2:3] == expected[2:3]
                @test c.stage_error ≈ cr.stage_error rtol = 1.0e-12
                @test c.tmp_embedded ≈ cr.tmp_embedded rtol = 1.0e-12
                @test c.Mtmp ≈ cr.Mtmp rtol = 1.0e-12
                @test c.Jk == stored
                @test i.stats.njacs == njacs && c.ncheckjac == ncheck
            end
            # A stage-system cache is not proof that the coefficient cache
            # is valid. Poisoned, invalid Jk must not be borrowed by the check.
            c.jac_valid = false
            c.stage_jac_valid = true
            c.Jk .= NaN
            njacs, ncheck = i.stats.njacs, c.ncheckjac
            tc = pade ? 0.12 : 0.1 + c.μ * 0.02
            distance = pade ? 0.02 : -c.μ * 0.02
            expected = ODTI.check_stage_defect!(reference, cr, 0.1, tc, distance, uprev, weights)
            observed = ODTI.check_stage_defect!(i, c, 0.1, tc, distance, uprev, weights)
            @test observed[1] ≈ expected[1] rtol = 1.0e-12
            @test c.Mtmp ≈ cr.Mtmp rtol = 1.0e-12
            @test i.stats.njacs == njacs + 1 && c.ncheckjac == ncheck + 1
            @test !c.jac_valid
        end
    end

    @testset "Discontinuities and reinit invalidate the borrowed derivatives" begin
        square!(du, u, p, t) = (du .= u .^ 2; nothing)
        p = ODEProblem(square!, [0.1], (0.0, 0.1))
        for make in (TaylorGaussA4, ImplicitTaylor2A4)
            i = init(p, make(); dt = 0.01, abstol = 1.0e-8, reltol = 1.0e-6)
            c = i.cache
            @test ODTI.implicit_taylor_stage!(i, c, 0.0, 0.01, [0.1], i.u)
            ncoef = c.ncoefjac
            SciMLBase.derivative_discontinuity!(i, true)
            @test ODTI.implicit_taylor_stage!(i, c, 0.01, 0.01, [0.2], i.u)
            @test c.ncoefjac > ncoef
            SciMLBase.reinit!(i; reset_dt = false)
            @test !c.jac_valid && !c.lu_valid && !c.filter_valid
        end
    end

    @testset "Time-varying nonnormal systems, scaling, and backward integration" begin
        # Exact solution u = scale * [cos(t), sin(t)]. Both the eigenvalues and
        # the eigenvectors of J(t) vary, so a good old contraction is no promise
        # that its linearization remains good on the next step.
        function moving!(du, u, p, t)
            scale = p
            a = 2 + sin(7t)
            x, y = u[1] - scale * cos(t), u[2] - scale * sin(t)
            du[1] = -a * x + 20cos(3t) * y - scale * sin(t)
            du[2] = -(4 + cos(5t)) * y + scale * cos(t)
            return nothing
        end
        for direction in (-1, 1), scale in (1.0e-6, 1.0, 1.0e6)
            t0, tf = direction == 1 ? (0.0, 0.5) : (0.5, 0.0)
            u0 = scale .* [cos(t0), sin(t0)]
            p = ODEProblem(moving!, u0, (t0, tf), scale)
            for make in (TaylorGaussA4, TaylorRadauLA5)
                sol = solve(p, make(); abstol = 1.0e-10scale, reltol = 1.0e-7)
                @test SciMLBase.successful_retcode(sol)
                @test maximum(
                    norm(u ./ scale - [cos(t), sin(t)], Inf)
                        for (u, t) in zip(sol.u, sol.t)
                ) < 2.0e-6
            end
        end
    end
end
