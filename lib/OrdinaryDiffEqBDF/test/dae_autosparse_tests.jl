using OrdinaryDiffEqBDF, SparseArrays, Test
import ADTypes, SciMLBase
using ADTypes: AutoFiniteDiff, AutoForwardDiff, AutoSparse
using SciMLBase: DAEFunction, DAEProblem, NoInit, isinplace, successful_retcode
using SparseConnectivityTracer: TracerSparsityDetector
using SparseMatrixColorings: GreedyColoringAlgorithm

const sparse_fd = AutoSparse(
    AutoForwardDiff(); sparsity_detector = TracerSparsityDetector(),
    coloring_algorithm = GreedyColoringAlgorithm()
)

# ∂F/∂du = diag(1, 1, 0) and ∂F/∂u are both strict subsets of the pattern of J.
function oscillator!(res, du, u, p, t)
    res[1] = du[1] - u[2]
    res[2] = du[2] + u[1]
    res[3] = u[3] - u[1]
    return nothing
end
oscillator(du, u, p, t) = [du[1] - u[2], du[2] + u[1], u[3] - u[1]]

# Heat equation with algebraic Dirichlet end nodes and a long-range du coupling,
# so ∂F/∂du has entries outside ∂F/∂u's tridiagonal band.
function heat!(res, du, u, p, t)
    n = length(u)
    res[1] = u[1] - sin(t)
    for i in 2:(n - 1)
        res[i] = du[i] - (u[i - 1] - 2u[i] + u[i + 1])
    end
    res[2] += du[n - 1] / 10
    res[n] = u[n]
    return nothing
end

@testset "DFBDF AutoSparse with differing ∂F/∂u and ∂F/∂du patterns" begin
    heat_n = 12
    problems = (
        DAEProblem(
            oscillator!, [0.0, -1.0, 0.0], [1.0, 0.0, 1.0], (0.0, 1.0);
            differential_vars = [true, true, false]
        ),
        DAEProblem{false}(
            oscillator, [0.0, -1.0, 0.0], [1.0, 0.0, 1.0], (0.0, 1.0);
            differential_vars = [true, true, false]
        ),
        DAEProblem(
            heat!, [1.0; zeros(heat_n - 1)], zeros(heat_n), (0.0, 1.0);
            differential_vars = [false; trues(heat_n - 2); false]
        ),
    )
    for prob in problems
        sol_dense = solve(prob, DFBDF(); abstol = 1.0e-10, reltol = 1.0e-10, initializealg = NoInit())
        sol_sparse = solve(
            prob, DFBDF(; autodiff = sparse_fd);
            abstol = 1.0e-10, reltol = 1.0e-10, initializealg = NoInit()
        )
        @test successful_retcode(sol_sparse)
        @test sol_sparse.u[end] ≈ sol_dense.u[end] rtol = 1.0e-8
        @test sol_sparse.stats.njacs > 0
        if isinplace(prob)
            integrator = init(prob, DFBDF(; autodiff = sparse_fd); initializealg = NoInit())
            @test integrator.cache.nlsolver.cache.J isa SparseMatrixCSC
            @test integrator.cache.nlsolver.cache.dae_jacobians.J_du isa SparseMatrixCSC
        end
    end
    @test solve(
        problems[1], DFBDF(; autodiff = sparse_fd); abstol = 1.0e-10, reltol = 1.0e-10
    ).u[end][1] ≈ cos(1) rtol = 1.0e-6
end

function decay!(r, du, u, p, t)
    r[1] = du[1] + u[1]
    r[2] = u[2] - u[1]
    return nothing
end

@testset "Known sparsity pattern with stored zeros" begin
    # The pattern's entries are all stored zeros; they are still structural.
    P = sparse([1, 2, 2], [1, 1, 2], [0.0, 0.0, 0.0], 2, 2)
    backend = AutoSparse(AutoFiniteDiff(); sparsity_detector = ADTypes.KnownJacobianSparsityDetector(P), coloring_algorithm = GreedyColoringAlgorithm())
    prob = DAEProblem(decay!, [-1.0, -1.0], [1.0, 1.0], (0.0, 0.4); differential_vars = [true, false])
    sol = solve(prob, DFBDF(autodiff = backend); initializealg = SciMLBase.NoInit(), abstol = 1.0e-10, reltol = 1.0e-10)
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end] ≈ fill(exp(-0.4), 2) atol = 1.0e-7 rtol = 1.0e-7

    integrator = init(prob, DFBDF(autodiff = backend); initializealg = NoInit())
    nlcache = integrator.cache.nlsolver.cache
    pattern(A) = findnz(A)[1:2]
    @test pattern(nlcache.J) == pattern(P)
    @test pattern(nlcache.dae_jacobians.J_du) == pattern(P)
    @test pattern(nlcache.W) == pattern(P)
    solve!(integrator)
    # After a step, `J` holds ∂F/∂u and `J_du` holds ∂F/∂(du).
    @test Matrix(nlcache.J) ≈ [1.0 0.0; -1.0 1.0]
    @test Matrix(nlcache.dae_jacobians.J_du) ≈ [1.0 0.0; 0.0 0.0] atol = 1.0e-12

    # Storage that cannot hold the known pattern must fail rather than skip derivatives.
    f_short = DAEFunction(decay!; jac_prototype = sparse([1, 2], [1, 2], [1.0, 1.0]))
    prob_short = DAEProblem(f_short, [-1.0, -1.0], [1.0, 1.0], (0.0, 0.4); differential_vars = [true, false])
    @test_throws DimensionMismatch init(prob_short, DFBDF(autodiff = backend); initializealg = NoInit())
end
