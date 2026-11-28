using OrdinaryDiffEqBDF, LinearSolve, ADTypes, SciMLBase, Test
import OrdinaryDiffEqDifferentiation
using SciMLOperators: MatrixOperator
using LinearAlgebra: I, Diagonal, mul!
using ADTypes: AutoForwardDiff, AutoFiniteDiff

# Forward differences destabilize stiff DAEs through the matrix-free operator, so the
# finite-difference tests use the centered stencil.
const centered_fd = AutoFiniteDiff(; fdtype = Val(:central))
using SciMLBase: DAEProblem, successful_retcode

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/1673

function decay!(res, du, u, p, t)
    @. res = du + u
    return nothing
end

# Semi-explicit index-1 DAE: u[3] is algebraic, so ∂F/∂du is singular.
function semi_explicit!(res, du, u, p, t)
    res[1] = du[1] + u[1] - u[2]^2
    res[2] = du[2] + 2u[2] - u[3]
    res[3] = u[3] - sin(u[1])
    return nothing
end

@testset "DFBDF with a matrix-free linear solver on a DAEProblem" begin
    prob = DAEProblem(decay!, -ones(4), ones(4), (0.0, 1.0))
    sol = solve(prob, DFBDF(; linsolve = KrylovJL_GMRES()); abstol = 1.0e-8, reltol = 1.0e-8)
    @test successful_retcode(sol)
    @test sol.u[end] ≈ fill(exp(-1), 4) rtol = 1.0e-5
    # Matrix-free: no Jacobian is ever formed.
    @test sol.stats.njacs == 0

    u0 = [1.0, 0.5, sin(1.0)]
    du0 = [-1.0 + 0.25, -1.0 + sin(1.0), 0.0]
    prob = DAEProblem(
        semi_explicit!, du0, u0, (0.0, 2.0); differential_vars = [true, true, false]
    )
    opts = (; abstol = 1.0e-9, reltol = 1.0e-9)
    reference = solve(prob, DFBDF(); opts...)
    for autodiff in (AutoForwardDiff(), centered_fd)
        sol = solve(prob, DFBDF(; autodiff, linsolve = KrylovJL_GMRES()); opts...)
        @test successful_retcode(sol)
        @test sol.stats.njacs == 0
        @test sol.u[end] ≈ reference.u[end] rtol = 1.0e-6
    end

    # User `precs` replace the residual weights; the Krylov stopping test is relative to
    # the residual in whatever norm they induce, so it holds without the weights too.
    linsolve = KrylovJL_GMRES(; precs = Returns((I, I)))
    sol = solve(prob, DFBDF(; linsolve); opts...)
    @test successful_retcode(sol)
    @test sol.u[end] ≈ reference.u[end] rtol = 1.0e-6
end

# The Krylov stopping test must scale with the residual: a residual scaled by 1e-12 once
# satisfied a fixed `abstol` at `dz = 0` and returned `Success` with `u = u0`.
f!(r, d, u, p, t) = (@. r = 1.0e-12 * (d + u); nothing)
p = DAEProblem(f!, -ones(4), ones(4), (0.0, 1.0))
@testset "Residual scale must not change solution" begin
    for ad in (AutoForwardDiff(), centered_fd)
        sol = solve(p, DFBDF(; linsolve = KrylovJL_GMRES(), autodiff = ad); abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end] ≈ fill(exp(-1), 4) rtol = 1.0e-5
    end
end

# Nonlinear in `du`, with exact solution `[exp(-t), exp(-2t)]` (the second factor of the
# first residual is positive).
function manufactured!(r, d, u, p, t)
    r[1] = (d[1] + u[1]) * (1 + d[1]^2 + u[2]^2)
    r[2] = u[2] - u[1]^2
    return nothing
end

@testset "Matrix-free DAE Newton solve with a large cj" begin
    prob = DAEProblem(
        manufactured!, [-1.0, -2.0], [1.0, 1.0], (0.0, 1.0); differential_vars = [true, false]
    )
    for autodiff in (AutoForwardDiff(), centered_fd)
        sol = solve(
            prob, DFBDF(; autodiff, linsolve = KrylovJL_GMRES()); abstol = 1.0e-7, reltol = 1.0e-7
        )
        @test successful_retcode(sol)
        @test maximum(maximum(abs, u .- [exp(-t), exp(-2t)]) for (u, t) in zip(sol.u, sol.t)) <
            1.0e-6
    end
end

# The matrix-free Newton test measures the Krylov residual through `Pl` applied by `ldiv!`.
# Under `ldiv = false`, with Krylov.jl's own `M`, `W` or `atol`, or with CG or MINRES, which
# measure `sqrt(r' * (Pl \ r))`, the solve tests another norm, and `ldiv = false` returned
# `Success` with `u = u0`. IterativeSolvers.jl wrappers are rejected because
# `LinearSolve.update_tolerances!` throws on them.
function swapped!(res, du, u, p, t)
    res[1] = du[2] + 2u[1] - sin(t)
    res[2] = du[1] + 2u[2] - cos(t)
    return nothing
end

@testset "Krylov solvers on the matrix-free DAE path" begin
    prob = DAEProblem(swapped!, [0.0, 0.0], [0.0, 0.5], (0.0, 1.0); differential_vars = [true, true])
    opts = (; abstol = 1.0e-8, reltol = 1.0e-8)
    reference = solve(prob, DFBDF(); opts...)
    for linsolve in (
            KrylovJL_GMRES(; ldiv = false), KrylovJL_CG(),
            KrylovJL_MINRES(),
            # What `IterativeSolversJL_BICGSTAB()` builds once IterativeSolvers is loaded.
            LinearSolve.IterativeSolversJL(identity, 0, nothing, (), (;)),
            KrylovJL_GMRES(; W = 1.0e-20I), KrylovJL_GMRES(; M = 1.0e20I),
            KrylovJL_GMRES(; atol = 1.0e20), SimpleGMRES(),
        )
        @test_throws ArgumentError solve(prob, DFBDF(; linsolve); opts...)
    end
    for linsolve in (KrylovJL_GMRES(), KrylovJL_FGMRES(), KrylovJL_BICGSTAB())
        sol = solve(prob, DFBDF(; linsolve); opts...)
        @test successful_retcode(sol)
        @test sol.u[end] ≈ reference.u[end] rtol = 1.0e-5
    end
end

# One fixed step of DFBDF is implicit Euler: `(u₁ - u₀) / dt + u₁ = 0`. The Krylov tolerance
# is measured in the left preconditioner the solve uses, and a Krylov solver can report
# convergence at `dz = 0` when that preconditioner is badly scaled. Whether the step was
# accepted is read from `integrator.t`, since `integrator.sol.retcode` is only set when the
# solve is finalized.
@testset "Krylov preconditioners and the implicit Euler step" begin
    prob = DAEProblem(decay!, [-1.0], [1.0], (0.0, 0.1))
    function euler_step(precs, warm_start)
        linsolve = KrylovJL_GMRES(; precs, warm_start)
        integrator = init(
            prob, DFBDF(; linsolve); dt = 0.1, adaptive = false,
            initializealg = SciMLBase.NoInit()
        )
        step!(integrator)
        u₁ = only(integrator.u)
        return integrator.t == 0.1, (u₁ - 1) / 0.1 + u₁
    end

    # `precs` returns a scalar diagonal until it sees the actual initial derivative,
    # which only the system being solved carries.
    changing = (A, ctx) -> (Diagonal([iszero(only(ctx[1])) ? 1.0e-20 : 1.0]), I)
    for warm_start in (LinearSolve.WarmStart.Auto, LinearSolve.WarmStart.None)
        accepted, residual = euler_step(changing, warm_start)
        @test accepted
        @test abs(residual) < 1.0e-12

        # Accepting `u = u₀` is the failure; a rejected step is acceptable.
        accepted, residual = euler_step(Returns((Diagonal([1.0e20]), I)), warm_start)
        @test !accepted
    end
end

# The Newton correction here is 1e-10, so a finite-difference product with it rather than
# with a unit vector perturbs `u = 1` by less than a rounding error.
@testset "Finite-difference validation of a tiny Newton correction" begin
    algebraic!(res, du, u, p, t) = (res[1] = (u[1] - 1) - t; nothing)
    prob = DAEProblem(algebraic!, [0.0], [1.0], (0.0, 1.0e-10); differential_vars = [false])
    integrator = init(
        prob, DFBDF(; autodiff = centered_fd, linsolve = KrylovJL_GMRES());
        dt = 1.0e-10, adaptive = false, initializealg = SciMLBase.NoInit()
    )
    step!(integrator)
    @test integrator.t == 1.0e-10
    @test only(integrator.u) ≈ 1 + 1.0e-10 atol = 1.0e-14
end

# A vector `reltol` once skipped the scaled tolerance and the true-residual check, and a
# preconditioner scaled by `1e20` made the solve return `Success` at `u = u0`. The linear
# solve takes the smallest component as its scalar tolerance.
@testset "Vector reltol on the matrix-free DAE path" begin
    prob = DAEProblem(decay!, -ones(2), ones(2), (0.0, 1.0))
    opts = (;
        dt = 0.1, abstol = 1.0e-8, reltol = fill(1.0e-8, 2),
        initializealg = SciMLBase.NoInit(), maxiters = 200,
    )
    badly_scaled = KrylovJL_GMRES(;
        precs = Returns((Diagonal(fill(1.0e20, 2)), I)), warm_start = LinearSolve.WarmStart.None
    )
    sol = solve(prob, DFBDF(; linsolve = badly_scaled); opts...)
    @test !successful_retcode(sol) || isapprox(sol.u[end], fill(exp(-1), 2); atol = 1.0e-5)

    sol = solve(
        prob, DFBDF(; linsolve = KrylovJL_GMRES()); abstol = 1.0e-8, reltol = [1.0e-8, 1.0e-6]
    )
    @test successful_retcode(sol)
    @test sol.u[end] ≈ fill(exp(-1), 2) atol = 1.0e-5
end

# A zero component of a vector `reltol` once made both Krylov stopping tolerances zero. The
# linear tolerance is floored at working precision; the integration tolerances stay
# componentwise.
function scaled_decay!(r, du, u, p, t)
    @. r = du + (1:3) * u
    return nothing
end

@testset "Vector reltol with a zero component" begin
    prob = DAEProblem(scaled_decay!, -[1.0, 2.0, 3.0], ones(3), (0.0, 1.0))
    exact = exp.(-[1.0, 2.0, 3.0])
    for linsolve in (KrylovJL_GMRES(), KrylovJL_FGMRES(), KrylovJL_BICGSTAB()),
            autodiff in (AutoForwardDiff(), centered_fd)
        autodiff isa AutoFiniteDiff && linsolve isa LinearSolve.KrylovJL &&
            nameof(linsolve.KrylovAlg) == :bicgstab! && continue
        sol = solve(
            prob, DFBDF(; linsolve, autodiff); abstol = 1.0e-8, reltol = [0.0, 1.0e-7, 1.0e-7]
        )
        @test successful_retcode(sol)
        @test sol.u[end] ≈ exact atol = 1.0e-6
        @test sol.stats.njacs == 0
    end
end

# Forward differences destabilize the matrix-free operator on stiff DAEs, and BiCGStab fails
# with the centered stencil, so finite differences are limited to the centered stencil with
# GMRES and FGMRES (issue #1673's Robertson DAE).
function robertson!(r, du, u, p, t)
    r[1] = du[1] + 0.04u[1] - 1.0e4 * u[2] * u[3]
    r[2] = du[2] - 0.04u[1] + 1.0e4 * u[2] * u[3] + 3.0e7 * u[2]^2
    r[3] = u[1] + u[2] + u[3] - 1
    return nothing
end

@testset "Robertson DAE with a matrix-free Krylov solve" begin
    prob = DAEProblem(
        robertson!, [-0.04, 0.04, 0.0], [1.0, 0.0, 0.0], (0.0, 1.0e5);
        differential_vars = [true, true, false]
    )
    reference = solve(prob, DFBDF(); abstol = 1.0e-12, reltol = 1.0e-10)
    supported = (
        (KrylovJL_GMRES(), AutoForwardDiff()), (KrylovJL_FGMRES(), AutoForwardDiff()),
        (KrylovJL_BICGSTAB(), AutoForwardDiff()),
        (KrylovJL_GMRES(), centered_fd), (KrylovJL_FGMRES(), centered_fd),
    )
    for (linsolve, autodiff) in supported
        sol = solve(prob, DFBDF(; linsolve, autodiff); abstol = 1.0e-9, reltol = 1.0e-6)
        @test successful_retcode(sol)
        @test sol.t[end] == 1.0e5
        @test sol.stats.njacs == 0
        @test maximum(abs(sum(u) - 1) for u in sol.u) < 1.0e-5
        @test maximum(abs, sol.u[end] - reference.u[end]) < 1.0e-5
    end
    unsupported = (
        (KrylovJL_GMRES(), AutoFiniteDiff()), (KrylovJL_FGMRES(), AutoFiniteDiff()),
        (KrylovJL_BICGSTAB(), AutoFiniteDiff()), (KrylovJL_BICGSTAB(), centered_fd),
        (KrylovJL_GMRES(), AutoFiniteDiff(; fdtype = Val(:complex))),
    )
    for (linsolve, autodiff) in unsupported
        @test_throws ArgumentError init(prob, DFBDF(; linsolve, autodiff))
    end
end

# Rows scaled far below the others are invisible to a residual norm; a preconditioner that
# equilibrates them (here the exact row scaling) restores the solve, whatever the scaling.
@testset "Algebraic constraint scaling on the matrix-free DAE path" begin
    for scale in (1.0, 1.0e-12)
        prob = DAEProblem(
            (r, du, u, p, t) -> (r[1] = du[1] + u[1] - t; r[2] = scale * (u[2] - t); nothing),
            [0.0, 0.0], [0.0, 0.0], (0.0, 0.1); differential_vars = [true, false]
        )
        linsolve = KrylovJL_GMRES(;
            precs = Returns((Diagonal([1.0, scale]), I)), warm_start = LinearSolve.WarmStart.None
        )
        sol = solve(
            prob, DFBDF(; linsolve); initializealg = SciMLBase.NoInit(), dt = 0.1,
            abstol = 1.0e-8, reltol = 1.0e-8
        )
        @test successful_retcode(sol)
        @test maximum(abs(u[2] - t) for (u, t) in zip(sol.u, sol.t)) <= 1.0e-8
    end
end

# A fixed step has no error control, so its Krylov forcing tolerance is the requested
# `reltol`, not machine epsilon, which GMRES cannot reach on a two-iteration solve.
@testset "Fixed-step implicit Euler with two Krylov iterations" begin
    prob = DAEProblem(
        (r, du, u, p, t) -> (r[1] = du[1] + u[1] - t; r[2] = u[2] - t; nothing),
        [0.0, 0.0], [0.0, 0.0], (0.0, 0.1); differential_vars = [true, false]
    )
    cold = KrylovJL_GMRES(; precs = Returns((I, I)), warm_start = LinearSolve.WarmStart.None)
    for linsolve in (KrylovJL_GMRES(), cold)
        sol = solve(
            prob, DFBDF(; linsolve); initializealg = SciMLBase.NoInit(), dt = 0.1, adaptive = false
        )
        @test successful_retcode(sol)
        @test sol.t[end] == 0.1
        @test sol.u[end] ≈ [0.1 / 11, 0.1] rtol = 1.0e-10
    end
end

# A zero Newton correction for a nonzero right-hand side is never accepted, however large the
# first residual (its allowance, `reltol * ‖P b₁‖`, can exceed the whole of a later residual).
# At `h = 1e-14` the first residual is `1e14`, the algebraic one `0.25`, and the first Krylov
# solve takes no iterations; the cold retry then supplies the correction. Each linear solve,
# a retry included, is counted in `stats.nsolve`.
@testset "Zero correction and linear solve count" begin
    for h in (1.0e-14, 1.0e-6)
        prob = DAEProblem(
            (r, du, u, p, t) -> (r[1] = du[1] - (t / h) / h; r[2] = u[2] - u[1]^2 / 4; nothing),
            [0.0, 0.0], [0.0, 0.0], (0.0, h)
        )
        trace = IOBuffer()
        linsolve = KrylovJL_GMRES(;
            precs = Returns((I, I)), warm_start = LinearSolve.WarmStart.None,
            verbose = 1, iostream = trace
        )
        sol = solve(
            prob, DFBDF(; linsolve); initializealg = SciMLBase.NoInit(), dt = h,
            adaptive = false, reltol = 1.0e-3, abstol = 1.0e-8
        )
        @test successful_retcode(sol)
        @test sol.t[end] == h
        @test sol.u[end] ≈ [1.0, 0.25] atol = 1.0e-8
        @test abs(sol.u[end][2] - sol.u[end][1]^2 / 4) <= 1.0e-8
        entries = count("GMRES", String(take!(trace)))
        # One more than the Newton iterations: the cold retry ran.
        @test entries == sol.stats.nsolve
        @test entries > sol.stats.nnonliniter
    end
end

# Acceptance is judged against the current iterate only: the first residual (here `1e14`)
# allows the Krylov solve to stop early, but must not excuse a retry that left the algebraic
# constraint unresolved. With one GMRES iteration the cold retry cannot satisfy it and the
# step is rejected; with two it can.
@testset "Retry that does not resolve the constraint is rejected" begin
    h = 1.0e-14
    prob = DAEProblem(
        (r, du, u, p, t) -> (
            r[1] = du[1] - (t / h) / h + u[1]^2 / 64; r[2] = u[2] - u[1]^2 / 4; nothing
        ),
        [0.0, 0.0], [0.0, 0.0], (0.0, h)
    )
    solve_with(itmax) = solve(
        prob,
        DFBDF(;
            linsolve = KrylovJL_GMRES(;
                precs = Returns((I, I)), warm_start = LinearSolve.WarmStart.None, itmax
            )
        );
        initializealg = SciMLBase.NoInit(), dt = h, adaptive = false, reltol = 1.0e-3,
        abstol = 1.0e-8
    )
    sol = solve_with(1)
    @test !successful_retcode(sol)
    @test sol.t[end] == 0.0

    sol = solve_with(2)
    @test successful_retcode(sol)
    @test sol.t[end] == h
    @test sol.u[end] ≈ [1.0, 0.25] atol = 1.0e-8
    @test abs(sol.u[end][2] - sol.u[end][1]^2 / 4) <= 1.0e-8
end

# A rejected cold retry still charges both linear solves.
@testset "Rejected retry is counted in stats.nsolve" begin
    prob = DAEProblem(decay!, [-1.0], [1.0], (0.0, 0.1))
    trace = IOBuffer()
    linsolve = KrylovJL_GMRES(;
        precs = Returns((Diagonal([1.0e20]), I)), warm_start = LinearSolve.WarmStart.None,
        verbose = 1, iostream = trace
    )
    sol = solve(
        prob, DFBDF(; linsolve); initializealg = SciMLBase.NoInit(), dt = 0.1, adaptive = false
    )
    @test sol.retcode == SciMLBase.ReturnCode.ConvergenceFailure
    @test sol.t[end] == 0.0
    entries = count("GMRES", String(take!(trace)))
    @test entries == 2
    @test sol.stats.nsolve == 2
end

# The acceptance test measures the residual in the left preconditioner. With `F` scaled by
# `1f-24` and `Pl = 1f24`, `Pl \ b` underflowed to zero, the zero correction passed and the
# step was accepted at `u = u₀`; an infinite `Pl` did the same at any scale. Any accepted step
# must solve the implicit Euler equation, and a preconditioner that equilibrates `F` (`Pl` =
# its scale) or a unit `F` must be accepted. Krylov.jl's GMRES returns `x = 0` for an
# operator of absolutely tiny norm, so a tiny `F` with `Pl = 1` may be rejected.
@testset "Preconditioned residual scale, underflow and nonfinite values" begin
    function scaled_euler_step(T, scale, pl; autodiff = AutoForwardDiff())
        F! = (r, du, u, p, t) -> (@. r = scale * (du + u); nothing)
        prob = DAEProblem(F!, T[-1], T[1], (zero(T), T(0.1)))
        linsolve = KrylovJL_GMRES(; precs = Returns((Diagonal(T[pl]), I)))
        integrator = init(
            prob, DFBDF(; linsolve, autodiff); dt = T(0.1), adaptive = false,
            initializealg = SciMLBase.NoInit()
        )
        step!(integrator)
        u₁ = only(integrator.u)
        return integrator.t == T(0.1), abs((u₁ - 1) / T(0.1) + u₁)
    end
    for T in (Float32, Float64), autodiff in (AutoForwardDiff(), centered_fd)
        tol = 100 * sqrt(eps(T))
        tiny, huge = T == Float32 ? (1.0f-24, 1.0f30) : (1.0e-80, 1.0e200)
        must_accept = [(one(T), one(T)), (T(1.0e12), one(T)), (T(1.0e12), T(1.0e-12))]
        may_reject = [(T(Inf), T(Inf)), (one(T), T(NaN))]
        for s in (one(T), T(1.0e-12), T(1.0e12), tiny, huge)
            push!(must_accept, (s, s))
            append!(may_reject, [(s, one(T)), (s, inv(s)), (s, T(Inf)), (s, T(NaN))])
        end
        # The reported cases: `Pl \ b` underflows (Float32) or `Pl = 1e280` (Float64).
        push!(may_reject, T == Float32 ? (1.0f-24, 1.0f24) : (1.0e-80, 1.0e280))
        for (s, pl) in must_accept
            accepted, residual = scaled_euler_step(T, s, pl; autodiff)
            @test accepted
            @test residual <= tol
        end
        for (s, pl) in may_reject
            accepted, residual = scaled_euler_step(T, s, pl; autodiff)
            @test !accepted || residual <= tol
        end
    end

    # A nonfinite residual is never accepted.
    nonfinite!(r, du, u, p, t) = (@. r = du + u + ifelse(t > 0, NaN, 0.0); nothing)
    prob = DAEProblem(nonfinite!, [-1.0], [1.0], (0.0, 0.1))
    integrator = init(
        prob, DFBDF(; linsolve = KrylovJL_GMRES()); dt = 0.1, adaptive = false,
        initializealg = SciMLBase.NoInit()
    )
    step!(integrator)
    @test integrator.t == 0.0

    # A subnormal operator equilibrated by `Pl` is an exact unit system: it must be solved.
    algebraic!(r, du, u, p, t) = (r[1] = 1.0e-309 * (u[1] - t); nothing)
    prob = DAEProblem(algebraic!, [0.0], [0.0], (0.0, 1.0); differential_vars = [false])
    linsolve = KrylovJL_GMRES(; precs = Returns((Diagonal([1.0e-309]), I)))
    integrator = init(
        prob, DFBDF(; linsolve); dt = 1.0, adaptive = false, initializealg = SciMLBase.NoInit()
    )
    step!(integrator)
    @test integrator.t == 1.0
    @test only(integrator.u) ≈ 1 rtol = 1.0e-10

    # Rows scaled apart by `F` and equilibrated by `Pl` are each resolved.
    for (T, s₁, s₂) in ((Float64, 1.0e12, 1.0e-12), (Float64, 1.0e150, 1.0e-150), (Float32, 1.0f15, 1.0f-15))
        mixed! = (r, du, u, p, t) -> (r[1] = s₁ * (du[1] + u[1]); r[2] = s₂ * (du[2] + u[2]); nothing)
        prob = DAEProblem(mixed!, T[-1, -1], T[1, 1], (zero(T), T(0.1)))
        linsolve = KrylovJL_GMRES(; precs = Returns((Diagonal(T[s₁, s₂]), I)))
        integrator = init(
            prob, DFBDF(; linsolve); dt = T(0.1), adaptive = false,
            initializealg = SciMLBase.NoInit()
        )
        step!(integrator)
        @test integrator.t == T(0.1)
        @test integrator.u ≈ fill(1 / T(1.1), 2) rtol = 100 * sqrt(eps(T))
    end
end

# The acceptance test applies `Pl` to the residual and right-hand side as they are: scaling
# them by `‖b‖` first erased the second row here, and the correction `[1, 0]`, whose
# left-preconditioned relative residual is about 0.71, was accepted.
@testset "Acceptance test resolves rows the preconditioner equilibrates" begin
    integrator = init(
        DAEProblem(decay!, [-1.0], [1.0], (0.0, 0.1)), DFBDF(; linsolve = KrylovJL_GMRES());
        initializealg = SciMLBase.NoInit()
    )
    Pl = Diagonal([1.0e300, 1.0e-100])
    linsolve = init(
        LinearProblem(Pl, [1.0e300, 1.0e-100]), KrylovJL_GMRES(; precs = Returns((Pl, I)))
    )
    solve!(linsolve)
    converged(u) = (
        linsolve.u .= u;
        OrdinaryDiffEqDifferentiation.linear_residual_converged!(integrator, linsolve, zeros(2), zeros(2))
    )
    @test !converged([1.0, 0.0])
    @test converged([1.0, 1.0])
    @test !converged([0.0, 0.0])

    # Norms near `floatmax`: the roundoff allowance once summed them before scaling, overflowed
    # and rejected an exact correction.
    for T in (Float32, Float64)
        big = T(0.9) * floatmax(T)
        exact = init(LinearProblem(Diagonal([big]), [big]), KrylovJL_GMRES())
        solve!(exact)
        exact.u .= one(T)
        @test OrdinaryDiffEqDifferentiation.linear_residual_converged!(
            integrator, exact, zeros(T, 1), zeros(T, 1)
        )
    end
    # With the identity `Pl` the step is exact. The default weight preconditioner makes
    # `Pl \ b` overflow, so GMRES itself starts from an infinite residual: the step may be
    # rejected, but never accepted wrong.
    huge!(r, du, u, p, t) = (r[1] = 1.0e308 * (u[1] - t); nothing)
    prob = DAEProblem(huge!, [0.0], [0.0], (0.0, 1.0); differential_vars = [false])
    for (linsolve, must_accept) in
        ((KrylovJL_GMRES(; precs = Returns((I, I))), true), (KrylovJL_GMRES(), false))
        huge_step = init(
            prob, DFBDF(; linsolve); dt = 1.0, adaptive = false,
            initializealg = SciMLBase.NoInit()
        )
        step!(huge_step)
        accepted = huge_step.t == 1.0
        must_accept && @test accepted
        @test !accepted || isapprox(only(huge_step.u), 1; rtol = 1.0e-10)
    end
end

# The matrix-free Newton matrix is built for the initial state dimension; every public change
# of the state size refuses before changing any state, so the integrator can carry on.
@testset "State-size changes of a matrix-free $(nameof(alg)) integrator" for (alg, rtol) in (
        (DImplicitEuler, 1.0e-3), (DABDF2, 1.0e-3), (DFBDF, 1.0e-5), (DNordsieckBDF, 1.0e-5),
    )
    prob = DAEProblem(decay!, -ones(2), ones(2), (0.0, 1.0))
    integrator = init(
        prob, alg(; linsolve = KrylovJL_GMRES()); abstol = 1.0e-8, reltol = 1.0e-8,
        initializealg = SciMLBase.NoInit()
    )
    step!(integrator, 0.5, true)
    for change! in (
            i -> resize!(i, 3), i -> resize!(i, 1), i -> resize!(i, (3,)),
            i -> deleteat!(i, 2), i -> addat!(i, 3:3),
        )
        @test_throws ArgumentError change!(integrator)
        @test length(integrator.u) == 2
        @test length(integrator.du) == 2
    end
    solve!(integrator)
    @test successful_retcode(integrator.sol)
    @test integrator.u ≈ fill(exp(-1), 2) rtol = rtol
end

# The matrix-free path never forms the Newton matrix: callbacks that write one are rejected at
# init, and a supplied `jac` or `jvp` is not called while stepping.
@testset "DAEFunction Jacobian fields on the matrix-free path" begin
    linsolve = KrylovJL_GMRES()
    for kw in (
            (; Wfact_t = (W, u, p, γ, t) -> nothing), (; Wfact = (W, u, p, γ, t) -> nothing),
        )
        prob = DAEProblem(DAEFunction(decay!; kw...), [-1.0], [1.0], (0.0, 1.0))
        @test_throws ArgumentError init(prob, DFBDF(; linsolve))
    end
    prob = DAEProblem(
        DAEFunction(
            decay!; jac = (J, du, u, p, γ, t) -> nothing,
            jac_prototype = MatrixOperator(zeros(1, 1))
        ), [-1.0], [1.0],
        (0.0, 1.0)
    )
    @test_throws ArgumentError init(prob, DFBDF(; linsolve))
    @test_throws ArgumentError init(prob, DFBDF())

    calls = Ref(0)
    jac!(J, du, u, p, γ, t) = (calls[] += 1; J .= (1 + γ) .* one.(J); nothing)
    jvp!(Jv, v, du, u, p, γ, t) = (calls[] += 1; @. Jv = (1 + γ) * v; nothing)
    prob = DAEProblem(DAEFunction(decay!; jac = jac!, jvp = jvp!), [-1.0], [1.0], (0.0, 1.0))
    sol = solve(
        prob, DFBDF(; linsolve); abstol = 1.0e-8, reltol = 1.0e-8,
        initializealg = SciMLBase.NoInit()
    )
    @test successful_retcode(sol)
    @test only(sol.u[end]) ≈ exp(-1) rtol = 1.0e-5
    @test calls[] == 0
end

# `cj` is a scalar, so the matrix-free `W` takes the current one every stage rather than when
# the step size moves past the cutoff that schedules a `precs` refresh. It once kept the
# first step's `cj` on a 10% step-size change.
@testset "Matrix-free W takes the current BDF coefficient" begin
    for autodiff in (AutoForwardDiff(), centered_fd)
        prob = DAEProblem(decay!, [-1.0], [1.0], (0.0, 1.0))
        integrator = init(
            prob, DFBDF(; linsolve = KrylovJL_GMRES(), autodiff, max_order = Val(1));
            dt = 1.0e-3, initializealg = SciMLBase.NoInit()
        )
        rtol = autodiff isa AutoFiniteDiff ? 1.0e-8 : 1.0e-12
        for h in (1.0e-3, 1.1e-3)
            set_proposed_dt!(integrator, h)
            step!(integrator)
            @test integrator.t - integrator.tprev ≈ h
            W = integrator.cache.nlsolver.cache.W
            @test only(mul!([0.0], W, [1.0])) ≈ 1 + 1 / (integrator.t - integrator.tprev) rtol = rtol
        end
        @test integrator.stats.njacs == 0
    end
end

# A centered finite-difference product evaluates both perturbed points of each of its two
# pushforwards; `stats.nf` once charged one evaluation per pushforward.
@testset "stats.nf counts every residual evaluation of a matrix-free DAE step" begin
    for autodiff in (AutoForwardDiff(), centered_fd)
        calls = Ref(0)
        counted!(r, du, u, p, t) = (calls[] += 1; @. r = du + u; nothing)
        prob = DAEProblem(counted!, [-1.0], [1.0], (0.0, 1.0))
        integrator = init(
            prob, DFBDF(; linsolve = KrylovJL_GMRES(), autodiff); dt = 0.01,
            initializealg = SciMLBase.NoInit()
        )
        calls₀, nf₀ = calls[], integrator.stats.nf
        step!(integrator)
        step!(integrator)
        @test calls[] - calls₀ == integrator.stats.nf - nf₀
    end
end
