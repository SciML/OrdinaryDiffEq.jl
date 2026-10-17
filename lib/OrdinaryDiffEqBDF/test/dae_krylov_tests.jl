using OrdinaryDiffEqBDF, LinearSolve, SparseArrays, Test
using LinearAlgebra: I
using ADTypes: AutoForwardDiff, AutoFiniteDiff, AutoSparse, KnownJacobianSparsityDetector
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

    # User `precs` replace the residual weights, so the weighted Krylov `abstol` must not
    # be applied to the unweighted solve. The tolerance is looser because an unweighted
    # solve stops at LinearSolve's default `abstol`, as it does on the ODE path.
    opts = (; abstol = 1.0e-6, reltol = 1.0e-6)
    reference = solve(prob, DFBDF(); opts...)
    linsolve = KrylovJL_GMRES(; precs = Returns((I, I)))
    sol = solve(prob, DFBDF(; linsolve); opts...)
    @test successful_retcode(sol)
    @test sol.u[end] ≈ reference.u[end] rtol = 1.0e-4
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
