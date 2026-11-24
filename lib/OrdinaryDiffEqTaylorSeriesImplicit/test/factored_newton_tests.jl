using LinearAlgebra

@testset "Implicit Taylor factored Newton matrix" begin
    ODTI = OrdinaryDiffEqTaylorSeriesImplicit

    # D(z) = 1 + Σ_k d_k z^k with d_k from the normalized increment weights.
    denominator(P, Q) = [ODTI.pade_increment_coefficients(Val(P), Val(Q))[2][k + 1] / factorial(k) for k in 1:Q]

    @testset "Factors multiply to the Padé denominator" begin
        rng = Random.Xoshiro(1)
        for (P, Q) in ((1, 3), (2, 3), (3, 3), (2, 4), (3, 5))
            d = denominator(P, Q)
            roots = ODTI.denominator_roots(P, Q)
            @test length(roots) == Q
            for _ in 1:5
                z = complex(randn(rng), randn(rng))
                @test prod(1 - z / r for r in roots) ≈ 1 + sum(d[k] * z^k for k in 1:Q) rtol = 1.0e-10
            end
            # Matrix form: the factored solve inverts D(hJ).
            n = 6
            J = randn(rng, n, n) - 3I
            h = 0.7
            F = ODTI.build_factored_newton(P, Q, Float64, n)
            @test ODTI.factor_denominator!(F, J, h)
            b = randn(rng, n)
            D = I + sum(d[k] * (h * J)^k for k in 1:Q)
            @test ODTI.solve_denominator!(F, copy(b)) ≈ D \ b rtol = 1.0e-10
            @test ODTI.nsolves(F) == Q
            @test isfinite(F.rounding[]) && F.rounding[] >= eps()
        end
    end

    @testset "1-norm estimate of the inverse from the LU factors" begin
        rng = Random.Xoshiro(2)
        for T in (Float64, ComplexF64), n in (1, 5, 40)
            for _ in 1:5
                A = randn(rng, T, n, n) + 2I
                exact = opnorm(inv(A), 1)
                LU, piv = copy(A), zeros(Int, n)
                @test ODTI.lu_factor!(LU, piv)
                est = ODTI.inverse_norm_estimate!(LU, zeros(T, n), zeros(T, n))
                # Hager's estimate is attained by some vector: a lower bound, and close.
                @test est <= exact * (1 + 1.0e-10)
                @test est >= exact / 3
                b = randn(rng, T, n)
                @test ODTI.lu_solve!(LU, piv, copy(b)) ≈ A \ b
            end
        end
    end

    @testset "Option plumbing" begin
        @test TaylorRadauLA5().newton_matrix === :assembled
        @test_throws ArgumentError TaylorRadauLA5(newton_matrix = :something)
        f!(du, u, p, t) = (du .= -u; nothing)
        prob = ODEProblem(f!, [1.0, 2.0], (0.0, 1.0))
        cache(alg) = ODTI.newton_state(init(prob, alg).cache).factored
        # :assembled keeps an inert object of the same type (the cache type does not
        # depend on the option value): size 0, never enabled or active.
        F0 = cache(TaylorRadauLA5())
        @test F0 isa ODTI.FactoredNewton && !F0.enabled && !F0.active[] && isempty(F0.correction)
        @test cache(TaylorRadauLA5(newton_matrix = :auto)) isa ODTI.FactoredNewton
        @test cache(TaylorLobattoLA6(newton_matrix = :auto)) isa ODTI.FactoredNewton
        # Q ≤ 2 keeps the assembled matrix (see the header of TaylorSeries_factored_newton.jl).
        @test cache(TaylorGaussA4()) === nothing
        @test cache(TaylorGaussA4(newton_matrix = :factored)) === nothing
        F = cache(TaylorRadauLA5(newton_matrix = :auto))
        @test F.auto && !F.active[]
        @test cache(TaylorRadauLA5(newton_matrix = :factored)).active[]
    end

    # Linear stiff problem with eigenvalues from -1 to -1e7 in a random orthogonal basis: f'' = 0,
    # so the factored matrix is the exact Newton matrix, while every row of M carries the stiff
    # modes and even the row-scaling invariant cond(M, x) ≈ |d_Q| (h|λ|)^Q exceeds 1/eps once the
    # steps grow past the fast modes.
    rng = Random.Xoshiro(3)
    n = 6
    λ = -(10.0 .^ range(0, 7; length = n))
    V = Matrix(qr(randn(rng, n, n)).Q)
    A = V * Diagonal(λ) * V'
    linear!(du, u, p, t) = (mul!(du, A, u); nothing)
    u0 = ones(n)
    stiff = ODEProblem(linear!, u0, (0.0, 50.0))
    function run_flags(prob, alg; kw...)
        i = init(prob, alg; save_everystep = false, kw...)
        c = ODTI.newton_state(i.cache)
        switched = false
        while i.t < prob.tspan[2]
            step!(i)
            switched |= c.factored !== nothing && c.factored.active[]
            i.sol.retcode in (SciMLBase.ReturnCode.Default, SciMLBase.ReturnCode.Success) || break
        end
        return i, switched
    end

    @testset "Switches when the assembled matrix is too ill-conditioned" begin
        for make in (TaylorRadauLA5, TaylorLobattoLA4, TaylorLobattoLA6)
            ia, sa = run_flags(stiff, make(newton_matrix = :assembled, stages = :reduced); abstol = 1.0e-8, reltol = 1.0e-6)
            i, s = run_flags(stiff, make(newton_matrix = :auto); abstol = 1.0e-8, reltol = 1.0e-6)
            @test !sa
            @test s
            @test SciMLBase.successful_retcode(i.sol)
            @test i.stats.naccept <= ia.stats.naccept
            @test norm(i.u - exp(50A) * u0, Inf) < 1.0e-5
        end
    end

    @testset "Forced factored matrix solves the linear problem" begin
        for make in (TaylorRadauLA5, TaylorLobattoLA4)
            sol = solve(stiff, make(newton_matrix = :factored); abstol = 1.0e-8, reltol = 1.0e-6)
            @test SciMLBase.successful_retcode(sol)
            @test norm(sol.u[end] - exp(50A) * u0, Inf) < 1.0e-5
        end
    end

    @testset "No switch when the ill-conditioning is row scaling" begin
        # The same eigenvalues on the diagonal: cond(M) is as large, but cond(M, x) = 1 and LU
        # with partial pivoting solves M exactly up to componentwise rounding.
        diag!(du, u, p, t) = (du .= λ .* u; nothing)
        dprob = ODEProblem(diag!, u0, (0.0, 50.0))
        ia, _ = run_flags(dprob, TaylorRadauLA5(newton_matrix = :assembled, stages = :reduced); abstol = 1.0e-8, reltol = 1.0e-6)
        i, s = run_flags(dprob, TaylorRadauLA5(newton_matrix = :auto); abstol = 1.0e-8, reltol = 1.0e-6)
        @test !s
        @test i.stats.naccept == ia.stats.naccept
    end

    @testset "No switch while the assembled matrix is well conditioned" begin
        # Moderate stiffness: h|λ| stays far below eps^(-1/Q).
        mild = ODEProblem(linear!, u0, (0.0, 1.0e-6))
        _, s = run_flags(mild, TaylorRadauLA5(newton_matrix = :auto); abstol = 1.0e-10, reltol = 1.0e-8)
        @test !s
        _, s = run_flags(ODEProblemLibrary.prob_ode_hires, TaylorRadauLA5(newton_matrix = :auto); abstol = 1.0e-9, reltol = 1.0e-7)
        @test !s
    end

    @testset "reinit returns to the assembled matrix" begin
        i, s = run_flags(stiff, TaylorRadauLA5(newton_matrix = :auto); abstol = 1.0e-8, reltol = 1.0e-6)
        @test s
        F = ODTI.newton_state(i.cache).factored
        F.active[] = true
        SciMLBase.reinit!(i)
        @test !F.active[] && isnan(F.distance[]) && F.rounding[] == Inf
        solve!(i)
        @test SciMLBase.successful_retcode(i.sol)
        @test norm(i.u - exp(50A) * u0, Inf) < 1.0e-5
    end
end
