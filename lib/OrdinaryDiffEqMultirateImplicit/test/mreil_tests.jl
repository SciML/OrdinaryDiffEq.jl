using OrdinaryDiffEqMultirateImplicit, OrdinaryDiffEqMultirate, OrdinaryDiffEqLowOrderRK,
    DiffEqDevTools, Test, LinearAlgebra
using SparseArrays, LinearSolve, SciMLOperators
using ADTypes: AutoFiniteDiff

# Stiff fast subsystem with a closed-form solution:
#   u1' = -λ (u1 - u2)   (fast, stiff)
#   u2' = -u2            (slow)
function stiff_split_problem(λ, u0, tspan)
    ffast!(du, u, p, t) = (du[1] = -λ * (u[1] - u[2]); du[2] = 0.0)
    fslow!(du, u, p, t) = (du[1] = 0.0; du[2] = -u[2])
    function analytic(u0, p, t)
        c = λ * u0[2] / (λ - 1)
        return [(u0[1] - c) * exp(-λ * t) + c * exp(-t), u0[2] * exp(-t)]
    end
    return SplitODEProblem(
            ODEFunction(ffast!; analytic = analytic), fslow!, u0, tspan
        ), analytic
end

# Manufactured solution y(t) = [exp(-t), sin(t)] whose fast Jacobian
# ∂f1/∂u = [-2 0; u2 u1] is state dependent and whose slow part depends on `t`,
# so both approximations of the base method — the frozen Jacobian and the frozen
# slow rate — are actually exercised.
function manufactured_split_problem()
    ffast!(du, u, p, t) = (du[1] = -2 * u[1] + exp(-t); du[2] = u[1] * u[2])
    fslow!(du, u, p, t) = (du[1] = 0.0; du[2] = cos(t) - exp(-t) * sin(t))
    analytic(u0, p, t) = [exp(-t), sin(t)]
    return SplitODEProblem(
        ODEFunction(ffast!; analytic = analytic), fslow!, [1.0, 0.0], (0.0, 1.0)
    )
end

# Large stiff fast subsystem with a banded (cyclic tridiagonal) Jacobian.
function banded_split_problem(N, tspan)
    function ffast!(du, u, p, t)
        @inbounds for i in 1:N
            du[i] = -100.0 * u[i] + 10.0 * u[mod1(i + 1, N)]
        end
        return du
    end
    fslow!(du, u, p, t) = (du .= -0.1 .* u)
    function jac!(J, u, p, t)
        @inbounds for i in 1:N
            J[i, i] = -100.0
            J[i, mod1(i + 1, N)] += 10.0
        end
        return J
    end
    prototype = spzeros(N, N)
    @inbounds for i in 1:N
        prototype[i, i] = 1.0
        prototype[i, mod1(i + 1, N)] = 1.0
    end
    return ffast!, fslow!, jac!, prototype, ones(N), tspan
end

@testset "MREIL" begin
    @testset "Construction" begin
        @test MREIL().m == 4
        @test MREIL().order == 4
        @test MREIL().seq === :harmonic
        @test MREIL(m = 8, order = 3, seq = :romberg).m == 8
        @test MREIL(m = 8, order = 3, seq = :romberg).order == 3
        @test MREIL(m = 8, order = 3, seq = :romberg).seq === :romberg
        @test_throws ArgumentError MREIL(m = 0)
        @test_throws ArgumentError MREIL(order = 1)
    end

    @testset "Scalar out-of-place" begin
        prob = SplitODEProblem(
            (u, p, t) -> -0.9 * u, (u, p, t) -> -0.1 * u, 1.0, (0.0, 1.0)
        )

        sol = solve(prob, MREIL(m = 10, order = 4), dt = 0.1, adaptive = false)
        @test abs(sol.u[end] - exp(-1.0)) < 1.0e-6

        sol_a = solve(prob, MREIL(m = 10, order = 4), reltol = 1.0e-8, abstol = 1.0e-8)
        @test abs(sol_a.u[end] - exp(-1.0)) < 1.0e-6
    end

    @testset "Vector in-place" begin
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        u0 = [1.0, 2.0, 3.0]
        prob = SplitODEProblem(f1!, f2!, u0, (0.0, 1.0))

        sol = solve(prob, MREIL(m = 10, order = 4), dt = 0.1, adaptive = false)
        @test norm(sol.u[end] - u0 .* exp(-1.0)) < 1.0e-6

        sol_a = solve(prob, MREIL(m = 10, order = 4), reltol = 1.0e-8, abstol = 1.0e-8)
        @test norm(sol_a.u[end] - u0 .* exp(-1.0)) < 1.0e-6

        sol_fd = solve(
            prob, MREIL(m = 10, order = 4, autodiff = AutoFiniteDiff()),
            dt = 0.1, adaptive = false
        )
        @test norm(sol_fd.u[end] - u0 .* exp(-1.0)) < 1.0e-6
    end

    @testset "Order convergence" begin
        f1_analytic = (u0, p, t) -> u0 * exp(-t)
        f1_inplace = (du, u, p, t) -> (du .= -0.9 .* u)
        f2_inplace = (du, u, p, t) -> (du .= -0.1 .* u)

        f1_func = ODEFunction(f1_inplace; analytic = f1_analytic)
        prob = SplitODEProblem(f1_func, f2_inplace, [1.0], (0.0, 1.0))

        dts = 1 .// 2 .^ (6:-1:2)
        testTol = 0.3

        for target_order in [2, 3, 4]
            sim = test_convergence(
                dts, prob, MREIL(m = 10, order = target_order), adaptive = false
            )
            @test sim.𝒪est[:final] ≈ target_order atol = testTol
        end

        for target_order in [2, 3, 4]
            sim = test_convergence(
                dts, prob, MREIL(m = 10, order = target_order, seq = :romberg),
                adaptive = false
            )
            @test sim.𝒪est[:final] ≈ target_order atol = testTol
        end
    end

    @testset "Order convergence, state- and time-dependent problem" begin
        prob = manufactured_split_problem()

        # Finer ladder than the linear problem above: order 4 is outside the
        # asymptotic regime at 1//2^2 here.
        dts = 1 .// 2 .^ (7:-1:4)
        testTol = 0.3

        for seq in (:harmonic, :romberg), target_order in [2, 3, 4]
            sim = test_convergence(
                dts, prob, MREIL(m = 4, order = target_order, seq = seq),
                adaptive = false
            )
            @test sim.𝒪est[:final] ≈ target_order atol = testTol
        end
    end

    @testset "Dense output" begin
        f1_func = ODEFunction(
            (du, u, p, t) -> (du .= -0.9 .* u); analytic = (u0, p, t) -> u0 * exp(-t)
        )
        f2_inplace = (du, u, p, t) -> (du .= -0.1 .* u)
        prob = SplitODEProblem(f1_func, f2_inplace, [1.0], (0.0, 1.0))

        sol = solve(prob, MREIL(m = 10, order = 4), reltol = 1.0e-10, abstol = 1.0e-10)
        @test sol.dense
        step_err = maximum(abs(sol.u[i][1] - exp(-sol.t[i])) for i in eachindex(sol.t))
        ts = range(0.0, 1.0, length = 997)
        dense_err = maximum(abs(sol(t)[1] - exp(-t)) for t in ts)
        @test step_err < 1.0e-10
        # Without k[1]/k[2] the Hermite interpolant is built from stale derivatives
        # and this is ~5e-3, eight orders above the error at the step points.
        @test dense_err < 1.0e-7

        # …and the interpolant converges at the method's order, not at order 1.
        errs = map((1 // 16, 1 // 32)) do dt
            s = solve(prob, MREIL(m = 4, order = 4), dt = dt, adaptive = false)
            maximum(abs(s(t)[1] - exp(-t)) for t in range(0.0, 1.0, length = 1001))
        end
        @test log2(errs[1] / errs[2]) > 3.5
    end

    @testset "Stiff fast subsystem" begin
        u0 = [1.0, 1.0]
        prob, analytic = stiff_split_problem(500.0, u0, (0.0, 1.0))
        exact = analytic(u0, nothing, 1.0)

        # λ * h_fast = 500 * 0.05 / 4 = 6.25, far outside explicit Euler's
        # stability region: MREEF diverges here while MREIL stays accurate.
        sol_il = solve(prob, MREIL(m = 4, order = 4), dt = 0.05, adaptive = false)
        sol_ef = solve(prob, MREEF(m = 4, order = 4), dt = 0.05, adaptive = false)
        @test norm(sol_il.u[end] - exact) < 1.0e-5
        @test norm(sol_ef.u[end] - exact) > 1.0e3

        # Adaptive: MREEF's step size is stability limited and grows with λ,
        # MREIL's is accuracy limited and stays put.
        prob4, analytic4 = stiff_split_problem(1.0e4, u0, (0.0, 1.0))
        exact4 = analytic4(u0, nothing, 1.0)
        sol_il4 = solve(prob4, MREIL(m = 4, order = 4), reltol = 1.0e-6, abstol = 1.0e-6)
        sol_ef4 = solve(prob4, MREEF(m = 4, order = 4), reltol = 1.0e-6, abstol = 1.0e-6)
        @test norm(sol_il4.u[end] - exact4) < 1.0e-5
        @test sol_il4.stats.naccept < sol_ef4.stats.naccept / 10
    end

    @testset "Two-way coupled stiff fast subsystem" begin
        # The slow rate reads the fast state here, unlike `stiff_split_problem`, so
        # the frozen slow rate is a real error source. Assert stability and step
        # count rather than a tight error bound: at fixed `dt` the frozen slow rate
        # puts an O(1/λ) floor under the error that extrapolation cannot remove.
        λ = 1.0e4
        ffast!(du, u, p, t) = (du[1] = -λ * (u[1] - u[2]); du[2] = 0.0)
        fslow!(du, u, p, t) = (du[1] = 0.0; du[2] = u[1] - u[2])
        prob = SplitODEProblem(ffast!, fslow!, [1.0, 0.5], (0.0, 1.0))

        ref = solve(prob, MREIL(m = 4, order = 4), reltol = 1.0e-12, abstol = 1.0e-12)
        sol_il = solve(prob, MREIL(m = 4, order = 4), reltol = 1.0e-8, abstol = 1.0e-8)
        sol_ef = solve(prob, MREEF(m = 4, order = 4), reltol = 1.0e-8, abstol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol_il)
        @test norm(sol_il.u[end] - ref.u[end]) < 1.0e-6
        @test sol_il.stats.naccept < sol_ef.stats.naccept / 10

        # Fixed dt: halving `dt` does not help once the fast transient is
        # unresolved, and the residual error tracks 1/λ.
        errs = map((0.05, 0.025)) do dt
            s = solve(prob, MREIL(m = 4, order = 4), dt = dt, adaptive = false)
            norm(s.u[end] - ref.u[end])
        end
        @test all(e -> e < 10 / λ, errs)
        @test errs[2] > errs[1] / 4
    end

    @testset "Nonlinear coupled system" begin
        function ff!(du, u, p, t)
            du[1] = 0.0
            du[2] = 20.0 * (2.0 * u[1] - u[1]^2 * u[2])
        end
        function fs!(du, u, p, t)
            du[1] = 1.0 + u[1]^2 * u[2] - 3.0 * u[1]
            du[2] = 0.0
        end
        prob = SplitODEProblem(ff!, fs!, [1.5, 3.0], (0.0, 0.5))
        ref = solve(prob, SplitEuler(), dt = 1.0e-6, adaptive = false)

        sol = solve(prob, MREIL(m = 20, order = 4), dt = 0.01, adaptive = false)
        @test norm(sol.u[end] - ref.u[end]) < 1.0e-3

        sol_a = solve(prob, MREIL(m = 20, order = 4), reltol = 1.0e-6, abstol = 1.0e-6)
        @test norm(sol_a.u[end] - ref.u[end]) < 1.0e-3
    end

    @testset "Stats tracking" begin
        # For MREIL(m, order, seq = :harmonic): ns[j] = j, so a step takes
        # sum(1:order) macro intervals, m * sum(1:order) linear solves, `order`
        # W factorizations (one per extrapolation column) and a single Jacobian.
        # Each step adds one f1 and one f2 evaluation at (u, t + dt) for k[2].
        # The Jacobian is supplied so that it costs no f1 evaluations of its own.
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        jac!(J, u, p, t) = (J[1, 1] = -0.9; J)
        prob = SplitODEProblem(
            ODEFunction(f1!; jac = jac!, jac_prototype = zeros(1, 1)), f2!,
            [1.0], (0.0, 0.5)
        )
        m_val, order = 5, 3
        sol = solve(
            prob, MREIL(m = m_val, order = order), dt = 0.1, adaptive = false
        )
        nsteps = sol.stats.naccept

        @test sol.stats.njacs == nsteps
        @test sol.stats.nw == order * nsteps
        @test sol.stats.nsolve == m_val * sum(1:order) * nsteps
        @test sol.stats.nf == (m_val * sum(1:order) + 1) * nsteps + 1
        @test sol.stats.nf2 == (sum(1:order) + 1) * nsteps + 1
    end

    @testset "Jacobian reuse across rejected steps" begin
        # J is taken at (uprev, t), which a rejected step does not change, so a
        # retry must reuse it. W still has to be rebuilt every column: dt changed.
        ffast!(du, u, p, t) = (
            du[1] = -0.04u[1]; du[2] = 0.04u[1] - 3.0e7 * u[2]^2; du[3] = 3.0e7 * u[2]^2
        )
        fslow!(du, u, p, t) = (
            du[1] = 1.0e4 * u[2] * u[3]; du[2] = -1.0e4 * u[2] * u[3]; du[3] = 0.0
        )
        prob = SplitODEProblem(ffast!, fslow!, [1.0, 0.0, 0.0], (0.0, 1.0))
        order = 4
        sol = solve(prob, MREIL(m = 4, order = order), reltol = 1.0e-6, abstol = 1.0e-8)

        @test sol.stats.nreject > 0
        @test sol.stats.njacs == sol.stats.naccept
        @test sol.stats.nw == order * (sol.stats.naccept + sol.stats.nreject)
    end

    @testset "Cached Jacobian is dropped whenever user code could have run" begin
        # Reuse is gated on (`iter`, `success_iter`), which singles out a retry of a
        # rejected attempt — the one case where no user code can have run, because
        # `handle_callbacks!` is only reached on the accept branch. Everything below
        # either follows an accepted step (so the next `loopheader!` advances
        # `success_iter`) or goes through `reinit!` (which zeroes both counters, while
        # a stamped `J_iter` is always ≥ 1 once a step has run). Note `reinit!` skips
        # `initialize!` under `reinit_cache = false`, so clearing a flag there would
        # not have covered that path, and keying on (uprev, t) would not either:
        # `reinit!` restores both.
        f1!(du, u, p, t) = (du .= p[1] .* u)
        f2!(du, u, p, t) = (du .= 0.0)
        mkinteg(; kw...) = init(
            SplitODEProblem(f1!, f2!, [1.0], (0.0, 1.0), [-1.0]),
            MREIL(m = 4, order = 4); dt = 1.0, adaptive = false,
            save_everystep = false, kw...
        )

        # A Jacobian taken at p = -1 and then used against p = -100 does not merely
        # lose accuracy, it diverges to ~1e12 instead of decaying to ~0, so `u` is an
        # independent witness alongside the exact `njacs` count.
        for (name, mutate!) in (
                ("p mutated between step! calls", integ -> nothing),
                ("reinit!", integ -> reinit!(integ, [1.0])),
                (
                    "reinit!(; reinit_cache = false)",
                    integ -> reinit!(integ, [1.0]; reinit_cache = false),
                ),
                (
                    "reinit!(; reinit_cache = false, erase_sol = false)",
                    integ -> reinit!(integ, [1.0]; reinit_cache = false, erase_sol = false),
                ),
                ("set_u!", integ -> set_u!(integ, [1.0])),
                ("set_t!", integ -> set_t!(integ, 0.0)),
                ("set_proposed_dt!", integ -> set_proposed_dt!(integ, 1.0)),
                ("add_tstop!", integ -> add_tstop!(integ, 5.0)),
                ("u_modified!", integ -> (set_u!(integ, [1.0]); u_modified!(integ, true))),
            )
            integ = mkinteg()
            step!(integ)
            integ.p[1] = -100.0
            mutate!(integ)
            njacs = integ.stats.njacs
            step!(integ)
            @test integ.stats.njacs == njacs + 1
            @test abs(integ.u[1]) < 1.0e-3
        end

        # A callback fires on the accept branch, so the attempt that follows it always
        # advances `success_iter`.
        fired = Ref(false)
        cb = DiscreteCallback(
            (u, t, integ) -> !fired[],
            integ -> (fired[] = true; integ.p[1] = -100.0)
        )
        integ = mkinteg(callback = cb, tstops = [0.5])
        step!(integ)
        njacs = integ.stats.njacs
        step!(integ)
        @test integ.stats.njacs == njacs + 1
    end

    @testset "SplitFunction wrapper" begin
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        ff = SplitFunction(f1!, f2!)
        prob = ODEProblem(ff, [1.0, 2.0], (0.0, 1.0))
        sol = solve(prob, MREIL(m = 10, order = 4), dt = 0.1, adaptive = false)
        @test norm(sol.u[end] - [1.0, 2.0] .* exp(-1.0)) < 1.0e-6
    end

    @testset "User-supplied fast Jacobian" begin
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        jac!(J, u, p, t) = (fill!(J, 0.0); J[diagind(J)] .= -0.9; J)
        prob = SplitODEProblem(
            ODEFunction(f1!; jac = jac!, jac_prototype = zeros(2, 2)), f2!,
            [1.0, 2.0], (0.0, 1.0)
        )
        sol = solve(prob, MREIL(m = 10, order = 4), dt = 0.1, adaptive = false)
        @test norm(sol.u[end] - [1.0, 2.0] .* exp(-1.0)) < 1.0e-6
    end

    @testset "Sparse fast Jacobian prototype" begin
        N = 64
        ffast!, fslow!, jac!, prototype, u0, tspan = banded_split_problem(N, (0.0, 0.05))
        prob_dense = SplitODEProblem(ffast!, fslow!, u0, tspan)
        prob_sparse = SplitODEProblem(
            ODEFunction(ffast!; jac = jac!, jac_prototype = copy(prototype)),
            fslow!, u0, tspan
        )

        sol_dense = solve(prob_dense, MREIL(m = 4, order = 4), reltol = 1.0e-10, abstol = 1.0e-10)
        sol_sparse = solve(prob_sparse, MREIL(m = 4, order = 4), reltol = 1.0e-10, abstol = 1.0e-10)
        @test norm(sol_sparse.u[end] - sol_dense.u[end]) / norm(sol_dense.u[end]) < 1.0e-10

        # `jac!` only touches stored entries, so the prototype's structure has to be
        # restored before every call (SciML/OrdinaryDiffEq.jl#2653): a second solve
        # with the same prototype must give bit-identical results.
        sol_again = solve(prob_sparse, MREIL(m = 4, order = 4), reltol = 1.0e-10, abstol = 1.0e-10)
        @test sol_again.u[end] == sol_sparse.u[end]
    end

    @testset "Jacobian prototype on the outer SplitFunction" begin
        # `nlsolve_f` resolves a SplitFunction to `f1`, whose `jac_prototype` is
        # `nothing` when the user set it on the outer function instead. The cache
        # merges the two before building `J`, so the step has to read that merged
        # prototype back off the cache; re-deriving it from `f1` skips the structure
        # reset above and `jac!`, which only accumulates into part of the pattern,
        # then grows without bound.
        N = 64
        ffast!, fslow!, jac!, prototype, u0, tspan = banded_split_problem(N, (0.0, 0.05))
        prob_outer = ODEProblem(
            SplitFunction(
                ODEFunction(ffast!; jac = jac!), ODEFunction(fslow!);
                jac_prototype = copy(prototype)
            ), u0, tspan
        )
        prob_inner = SplitODEProblem(
            ODEFunction(ffast!; jac = jac!, jac_prototype = copy(prototype)),
            fslow!, u0, tspan
        )

        # This Jacobian is constant, so every step has to reproduce the first one.
        integ = init(prob_outer, MREIL(m = 4, order = 4), dt = 0.005, adaptive = false)
        step!(integ)
        first_J = copy(integ.cache.J)
        for _ in 1:4
            step!(integ)
        end
        @test integ.cache.J == first_J

        # Where the prototype was written must not change the answer at all.
        alg = MREIL(m = 4, order = 4)
        @test solve(prob_outer, alg, dt = 0.005, adaptive = false).u[end] ==
            solve(prob_inner, alg, dt = 0.005, adaptive = false).u[end]
    end

    @testset "User Jacobian on the outer SplitFunction" begin
        # A `SplitFunction`'s `jac` is `df1/du` — SciMLBase documents the derivatives as
        # defined on the `f1` portion only, and the constructor defaults the field from
        # `f1` — so an outer `jac` is the *fast* Jacobian and must be used exactly as
        # the same function written on `f1` would be.
        A1 = [-20.0 9.0; 9.0 -15.0]
        A2 = [-3.0 1.5; 1.5 -2.5]
        ffast!(du, u, p, t) = mul!(du, A1, u)
        fslow!(du, u, p, t) = mul!(du, A2, u)
        ffast(u, p, t) = A1 * u
        fslow(u, p, t) = A2 * u
        jac!(J, u, p, t) = (J .= A1; J)
        jac(u, p, t) = copy(A1)
        u0 = [1.0, 2.0]
        tspan = (0.0, 0.5)
        alg = MREIL(m = 4, order = 4)

        on_outer = ODEProblem(
            SplitFunction(ODEFunction(ffast!), ODEFunction(fslow!); jac = jac!), u0, tspan
        )
        on_f1 = SplitODEProblem(ODEFunction(ffast!; jac = jac!), fslow!, u0, tspan)
        integ = init(on_outer, alg, dt = 0.01, adaptive = false)
        step!(integ)
        @test integ.cache.J == A1
        @test solve(on_outer, alg, dt = 0.01, adaptive = false).u[end] ==
            solve(on_f1, alg, dt = 0.01, adaptive = false).u[end]

        on_outer_oop = ODEProblem(
            SplitFunction(ODEFunction{false}(ffast), ODEFunction{false}(fslow); jac = jac),
            u0, tspan
        )
        on_f1_oop = SplitODEProblem(
            ODEFunction{false}(ffast; jac = jac), ODEFunction{false}(fslow), u0, tspan
        )
        @test solve(on_outer_oop, alg, dt = 0.01, adaptive = false).u[end] ==
            solve(on_f1_oop, alg, dt = 0.01, adaptive = false).u[end]

        # …and it is reached through the same resolved description `J` was built from.
        # Pairing the outer `jac` with the outer `jac_prototype` instead skips the
        # sparse structure reset, and a `jac` that only accumulates into stored entries
        # then grows every step (SciML/OrdinaryDiffEq.jl#2653).
        N = 8
        ffast2!, fslow2!, jac2!, prototype, u0b, tspan2 = banded_split_problem(N, (0.0, 0.05))
        split_sources = ODEProblem(
            SplitFunction(
                ODEFunction(ffast2!; jac_prototype = copy(prototype)), ODEFunction(fslow2!);
                jac = jac2!, jac_prototype = nothing
            ), u0b, tspan2
        )
        integ2 = init(split_sources, alg, dt = 0.005, adaptive = false)
        step!(integ2)
        first_J = copy(integ2.cache.J)
        for _ in 1:4
            step!(integ2)
        end
        @test integ2.cache.J == first_J
        @test first_J[1, 2] == 10.0
    end

    @testset "All-zero dense Jacobian prototype is rejected" begin
        # A dense prototype carries its sparsity pattern in its values, so an all-zero
        # one declares an empty pattern; AD then leaves J identically zero and the
        # linearly implicit substep silently becomes an explicit Euler substep. Since
        # `zeros(n, n)` is a natural way to ask for a dense n-by-n buffer, refuse it.
        ffast!(du, u, p, t) = (du .= -1000.0 .* u)
        fslow!(du, u, p, t) = (du .= -0.1 .* u)
        ffast(u, p, t) = -1000.0 .* u
        fslow(u, p, t) = -0.1 .* u
        u0 = [1.0, 2.0]
        tspan = (0.0, 1.0)
        exact = u0 .* exp(-1000.1)

        prob = SplitODEProblem(
            ODEFunction(ffast!; jac_prototype = zeros(2, 2)), fslow!, u0, tspan
        )
        @test_throws ArgumentError solve(prob, MREIL(), dt = 0.01, adaptive = false)

        prob_oop = SplitODEProblem(
            ODEFunction{false}(ffast; jac_prototype = zeros(2, 2)), fslow, u0, tspan
        )
        @test_throws ArgumentError solve(prob_oop, MREIL(), dt = 0.01, adaptive = false)

        # `SplitFunction` copies `sparsity` from `f1`, so reaching the same state from
        # the outer function takes writing `sparsity` there explicitly.
        prob_outer = ODEProblem(
            SplitFunction(
                ODEFunction(ffast!), ODEFunction(fslow!);
                jac_prototype = zeros(2, 2), sparsity = zeros(2, 2)
            ), u0, tspan
        )
        @test_throws ArgumentError solve(prob_outer, MREIL(), dt = 0.01, adaptive = false)

        # A user `jac` writes the prototype itself and never reads the pattern off it,
        # so an all-zero buffer is the normal thing to pass there.
        jac!(J, u, p, t) = (fill!(J, 0.0); J[diagind(J)] .= -1000.0; J)
        prob_jac = SplitODEProblem(
            ODEFunction(ffast!; jac = jac!, jac_prototype = zeros(2, 2)), fslow!, u0, tspan
        )
        @test norm(
            solve(prob_jac, MREIL(), dt = 0.01, adaptive = false).u[end] - exact
        ) < 1.0e-8

        # A sparse prototype states its pattern through its stored structure, so
        # all-zero stored values still name a pattern and must keep working.
        prob_sparse = SplitODEProblem(
            ODEFunction(ffast!; jac_prototype = SparseMatrixCSC(2, 2, [1, 2, 3], [1, 2], zeros(2))),
            fslow!, u0, tspan
        )
        @test norm(
            solve(prob_sparse, MREIL(), dt = 0.01, adaptive = false).u[end] - exact
        ) < 1.0e-8
    end

    @testset "Matrix-free linear solver" begin
        N = 32
        A = 100.0 .* Matrix(Tridiagonal(fill(1.0, N - 1), fill(-2.0, N), fill(1.0, N - 1)))
        f1!(du, u, p, t) = mul!(du, A, u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        u0 = ones(N)
        prob = SplitODEProblem(f1!, f2!, u0, (0.0, 0.1))
        exact = exp((A - 0.1I) * 0.1) * u0

        sol_lu = solve(
            prob, MREIL(m = 4, order = 4, linsolve = LUFactorization()),
            reltol = 1.0e-10, abstol = 1.0e-10
        )
        @test norm(sol_lu.u[end] - exact) / norm(exact) < 1.0e-8

        for alg in (
                MREIL(m = 4, order = 4, linsolve = KrylovJL_GMRES()),
                MREIL(m = 4, order = 4, linsolve = KrylovJL_GMRES(), concrete_jac = true),
            )
            sol = solve(prob, alg, reltol = 1.0e-10, abstol = 1.0e-10)
            @test SciMLBase.successful_retcode(sol)
            @test norm(sol.u[end] - exact) / norm(exact) < 1.0e-8
            @test norm(sol.u[end] - sol_lu.u[end]) / norm(exact) < 1.0e-8
        end
    end

    @testset "Semilinear (operator) fast component" begin
        N = 32
        A = 100.0 .* Matrix(Tridiagonal(fill(1.0, N - 1), fill(-2.0, N), fill(1.0, N - 1)))
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        u0 = ones(N)
        prob = SplitODEProblem(MatrixOperator(A), f2!, u0, (0.0, 0.1))
        exact = exp((A - 0.1I) * 0.1) * u0

        sol = solve(prob, MREIL(m = 4, order = 4), reltol = 1.0e-10, abstol = 1.0e-10)
        @test SciMLBase.successful_retcode(sol)
        @test norm(sol.u[end] - exact) / norm(exact) < 1.0e-8

        # W = J - M/h_fast is rebuilt per column even though J is constant, so a
        # fixed-dt run has to converge at the extrapolation order.
        errs = map((1 // 100, 1 // 200)) do dt
            s = solve(
                prob, MREIL(m = 4, order = 4, linsolve = LUFactorization()),
                dt = dt, adaptive = false
            )
            norm(s.u[end] - exact) / norm(exact)
        end
        @test log2(errs[1] / errs[2]) > 3.0
    end

    @testset "Linear stability function" begin
        # With f2 = 0 the fast substep collapses to (1 - hλ)Δ = hλu, i.e. plain
        # implicit Euler, so R(z) is the Aitken–Neville extrapolation of
        # (1 - z/N_j)^(-N_j) over the columns' substep counts N_j = n_j*m. The
        # docstring quotes numbers off this function; pin them, and pin the closed
        # form against the solver so the two cannot drift apart.
        function Rstab(z, m, order, seq)
            ns = seq === :harmonic ? [j for j in 1:order] : [1 << (j - 1) for j in 1:order]
            T = [(1 - z / (ns[j] * m))^(-ns[j] * m) for j in 1:order]
            for k in 1:(order - 1), j in order:-1:(k + 1)
                T[j] = T[j] + (T[j] - T[j - 1]) / (ns[j] / ns[j - k] - 1)
            end
            return T[order]
        end

        # u' = Au with A = [a -b; b a] realifies u' = (a + ib)u, so one step of size 1
        # from [1, 0] lands on (Re R(z), Im R(z)).
        function Rstab_solver(z, m, order)
            A = [real(z) -imag(z); imag(z) real(z)]
            prob = SplitODEProblem(
                (du, u, p, t) -> mul!(du, A, u), (du, u, p, t) -> (du .= 0.0),
                [1.0, 0.0], (0.0, 1.0)
            )
            integ = init(
                prob, MREIL(m = m, order = order), dt = 1.0,
                adaptive = false, save_everystep = false
            )
            step!(integ)
            return complex(integ.u[1], integ.u[2])
        end

        for order in 2:5, z in (-1.0 + 0.0im, -30.0 + 0.0im, 3.7im, -2.0 + 5.0im)
            @test Rstab(z, 4, order, :harmonic) ≈ Rstab_solver(z, 4, order) rtol = 1.0e-9
        end

        # R is analytic on Re z ≤ 0 (its poles sit at z = N_j > 0) and vanishes at
        # infinity, so by the maximum modulus principle the supremum over the closed
        # left half plane is attained on the imaginary axis. A-stability holds only at
        # order 2; above it |R| exceeds 1 there by the amounts the docstring quotes.
        supR(m, order, seq) = maximum(
            abs(Rstab(complex(0.0, y), m, order, seq))
                for y in range(0.0, 60.0, length = 60001)
        )
        @test supR(4, 2, :harmonic) <= 1 + 1.0e-12
        @test supR(4, 3, :harmonic) ≈ 1.000594 atol = 1.0e-6
        @test supR(4, 4, :harmonic) ≈ 1.001502 atol = 1.0e-6
        @test supR(4, 5, :harmonic) ≈ 1.002282 atol = 1.0e-6
        @test supR(8, 4, :harmonic) ≈ 1.000573 atol = 1.0e-6
        @test supR(4, 4, :romberg) ≈ 1.000808 atol = 1.0e-6

        # R(∞) = 0 at every order: each column's factor vanishes as |z| → ∞ and the
        # extrapolation weights sum to 1. That is the damping the linearly implicit
        # base method is wanted for, and it survives the loss of A-stability.
        for order in 2:5, seq in (:harmonic, :romberg)
            @test abs(Rstab(-1.0e6 + 0.0im, 4, order, seq)) < 1.0e-20
            @test abs(Rstab(1.0e6im, 4, order, seq)) < 1.0e-20
        end
    end

    @testset "Mass matrix" begin
        M = Diagonal([2.0, 3.0])
        exact = [exp(-1.0), 2 * exp(-1.0)]

        ffast!(du, u, p, t) = (du[1] = -2.0 * u[1]; du[2] = -3.0 * u[2])
        fslow!(du, u, p, t) = (du[1] = 0.0; du[2] = 0.0)
        prob = ODEProblem(
            SplitFunction(ODEFunction(ffast!), ODEFunction(fslow!); mass_matrix = M),
            [1.0, 2.0], (0.0, 1.0)
        )
        sol = solve(prob, MREIL(m = 4, order = 4), dt = 0.05, adaptive = false)
        @test norm(sol.u[end] - exact) < 1.0e-8

        ffast(u, p, t) = [-2.0 * u[1], -3.0 * u[2]]
        fslow(u, p, t) = [0.0, 0.0]
        prob_oop = ODEProblem(
            SplitFunction(
                ODEFunction{false}(ffast), ODEFunction{false}(fslow); mass_matrix = M
            ),
            [1.0, 2.0], (0.0, 1.0)
        )
        sol_oop = solve(prob_oop, MREIL(m = 4, order = 4), dt = 0.05, adaptive = false)
        @test norm(sol_oop.u[end] - exact) < 1.0e-8

        # Scalar state: `Diagonal`/`Matrix` mass matrices never reach the scalar W, so
        # only a `UniformScaling` exercises it. 2u' = -u, i.e. u(t) = exp(-t/2).
        prob_scalar = ODEProblem(
            SplitFunction(
                ODEFunction{false}((u, p, t) -> -0.9 * u),
                ODEFunction{false}((u, p, t) -> -0.1 * u); mass_matrix = 2.0I
            ),
            1.0, (0.0, 1.0)
        )
        sol_scalar = solve(prob_scalar, MREIL(m = 4, order = 4), dt = 0.05, adaptive = false)
        @test abs(sol_scalar.u[end] - exp(-0.5)) < 1.0e-8

        sol_scalar_a = solve(
            prob_scalar, MREIL(m = 4, order = 4), reltol = 1.0e-10, abstol = 1.0e-10
        )
        @test abs(sol_scalar_a.u[end] - exp(-0.5)) < 1.0e-8

        # An array mass matrix on a scalar state is rejected rather than dropped.
        prob_bad = ODEProblem(
            SplitFunction(
                ODEFunction{false}((u, p, t) -> -0.9 * u),
                ODEFunction{false}((u, p, t) -> -0.1 * u); mass_matrix = Diagonal([2.0])
            ),
            1.0, (0.0, 1.0)
        )
        @test_throws ArgumentError solve(
            prob_bad, MREIL(m = 4, order = 4), dt = 0.05, adaptive = false
        )
    end
end
