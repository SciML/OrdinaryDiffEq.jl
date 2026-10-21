using OrdinaryDiffEqBDF, LinearSolve, ADTypes, SciMLBase, SparseArrays, Test
using LinearAlgebra: I, mul!, norm
using ADTypes: AutoForwardDiff, AutoFiniteDiff, AutoSparse, KnownJacobianSparsityDetector
using SciMLBase: DAEProblem, DAEFunction, successful_retcode
using SciMLOperators: update_coefficients!
using OrdinaryDiffEqDifferentiation: dae_jvp_cache, set_dae_jvp_cj!

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
# first residual is positive). At `u = [1, 1]`, `du = [-1, -2]` its Newton matrix is
# `[3(cj + 1) 0; -2 1]`.
function manufactured!(r, d, u, p, t)
    r[1] = (d[1] + u[1]) * (1 + d[1]^2 + u[2]^2)
    r[2] = u[2] - u[1]^2
    return nothing
end

@testset "Matrix-free DAE Newton matrix at a large cj" begin
    u, du, v = [1.0, 1.0], [-1.0, -2.0], [1.0, 0.5]
    for autodiff in (AutoForwardDiff(), AutoFiniteDiff()), cj in (1.0, 1.0e6, 1.0e10)
        J = dae_jvp_cache(DAEFunction{true}(manufactured!), copy(u), nothing, 0.0; autodiff)
        update_coefficients!(J, copy(u), nothing, 0.0; du = copy(du))
        set_dae_jvp_cj!(J, cj)
        exact = [3(cj + 1) 0; -2 1] * v
        @test norm(mul!(zeros(2), J, v) - exact) <= 1.0e-6 * norm(exact)
    end

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

# A linear index-1 DAE whose residual rows are rotated against the unknowns, the way a
# ModelingToolkit discretization can order them: every diagonal entry of
# ∂F/∂u + cj * ∂F/∂(du) is either zero or small against cj.
function rotated_chain!(res, du, u, p, t)
    n = length(u)
    for i in 1:(n - 1)
        res[i + 1] = du[i] + 2u[i] - u[i + 1]
    end
    res[1] = u[n] - u[1] / 2
    return nothing
end

@testset "Matrix-free DAE Krylov matches residual rows to unknowns" begin
    n = 40
    u0 = [ones(n - 1); 0.5]
    du0 = [fill(-1.0, n - 2); -1.5; -0.5]
    prob = DAEProblem(
        rotated_chain!, du0, u0, (0.0, 1.0); differential_vars = [trues(n - 1); false]
    )
    opts = (; abstol = 1.0e-9, reltol = 1.0e-9)
    reference = solve(prob, DFBDF(); opts...)
    jac_pattern = sparse(
        [[i + 1 for i in 1:(n - 1)]; [i + 1 for i in 1:(n - 1)]; 1; 1],
        [1:(n - 1); 2:n; n; 1], true, n, n
    )
    sparse_fd = AutoSparse(
        AutoForwardDiff(); sparsity_detector = KnownJacobianSparsityDetector(jac_pattern)
    )
    iterations_per_solve = map((AutoForwardDiff(), sparse_fd)) do autodiff
        iterations = Ref(0)
        linsolve = KrylovJL_GMRES(; callback = _ -> (iterations[] += 1; false))
        sol = solve(prob, DFBDF(; autodiff, linsolve); opts...)
        @test successful_retcode(sol)
        @test sol.stats.njacs == 0
        @test sol.u[end] ≈ reference.u[end] rtol = 1.0e-7
        iterations[] / sol.stats.nsolve
    end
    # Without a sparsity pattern there is no matching, and GMRES needs about one iteration
    # per unknown; with one, the row permutation alone brings that down to a handful.
    @test iterations_per_solve[1] > n / 2
    @test iterations_per_solve[2] < 10

    # Rows already ordered like the unknowns match to the identity, which leaves the
    # default preconditioner, and so the solve, exactly as it is without a pattern.
    u0 = [1.0, 0.5, sin(1.0)]
    du0 = [-1.0 + 0.25, -1.0 + sin(1.0), 0.0]
    prob = DAEProblem(
        semi_explicit!, du0, u0, (0.0, 2.0); differential_vars = [true, true, false]
    )
    pattern = AutoSparse(
        AutoForwardDiff();
        sparsity_detector = KnownJacobianSparsityDetector(sparse(Bool[1 1 0; 0 1 1; 1 0 1]))
    )
    sols = map((AutoForwardDiff(), pattern)) do autodiff
        solve(prob, DFBDF(; autodiff, linsolve = KrylovJL_GMRES()); opts...)
    end
    @test sols[1].u == sols[2].u
end
