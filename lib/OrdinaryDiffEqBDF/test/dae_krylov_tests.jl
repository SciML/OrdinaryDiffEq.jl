using OrdinaryDiffEqBDF, LinearSolve, ADTypes, SciMLBase, Test
using LinearAlgebra: I, Diagonal
using ADTypes: AutoForwardDiff, AutoFiniteDiff
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
    for autodiff in (AutoForwardDiff(), AutoFiniteDiff())
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
    for ad in (AutoForwardDiff(), AutoFiniteDiff())
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
    for autodiff in (AutoForwardDiff(), AutoFiniteDiff())
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
            KrylovJL_GMRES(; atol = 1.0e20),
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
        prob, DFBDF(; autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES());
        dt = 1.0e-10, adaptive = false, initializealg = SciMLBase.NoInit()
    )
    step!(integrator)
    @test integrator.t == 1.0e-10
    @test only(integrator.u) ≈ 1 + 1.0e-10 atol = 1.0e-14
end
