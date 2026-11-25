using OrdinaryDiffEqBDF, OrdinaryDiffEqDifferentiation, LinearSolve, SparseArrays, Test
import ADTypes, SciMLBase
using ADTypes: AutoFiniteDiff, AutoForwardDiff, AutoSparse
using SciMLBase: DAEFunction, DAEProblem, NoInit, isinplace, successful_retcode
using SparseConnectivityTracer: TracerSparsityDetector, TracerLocalSparsityDetector
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
    # An explicit tracer detects the missing (2, 1) entry as well and must not be trusted
    # to fit the prototype either.
    for detector in (ADTypes.KnownJacobianSparsityDetector(P), TracerSparsityDetector())
        backend_short = AutoSparse(
            AutoForwardDiff(); sparsity_detector = detector,
            coloring_algorithm = GreedyColoringAlgorithm()
        )
        # A tracer is checked at the first Jacobian refresh, at the real evaluation point.
        @test_throws "does not hold" solve(
            prob_short, DFBDF(autodiff = backend_short); initializealg = NoInit()
        )
    end
end

# Rows 2 and 3 swap at t = 0.5, so the pattern of ∂F/∂u changes from I to the swap of
# columns 2 and 3 while the exact state exp.((-t, -2t, -3t)) is the same on both sides.
function switching!(r, du, u, p, t)
    r[1] = du[1] + u[1]
    a = u[2] - exp(-2t)
    b = u[3] - exp(-3t)
    r[2] = t < 0.5 ? a : u[3] - exp(-3t)
    r[3] = t < 0.5 ? b : u[2] - exp(-2t)
    return nothing
end

@testset "DAE AutoSparse with a pattern that changes in time" begin
    switching_prob = DAEProblem(
        switching!, [-1.0, -2.0, -3.0], [1.0, 1.0, 1.0], (0.0, 1.0);
        differential_vars = [true, false, false]
    )
    exact(t) = [exp(-t), exp(-2t), exp(-3t)]
    solve_kwargs = (; initializealg = NoInit(), dt = 0.01, tstops = [0.5])
    J_u_late = [1.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0]
    J_du = [1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    # The union of both patterns, which the caller supplies when it is known.
    union_pattern = sparse(
        [1, 2, 3, 2, 3], [1, 2, 3, 3, 2], ones(5), 3, 3
    )

    # An automatically detected pattern is frozen at the initial state, so the first
    # refresh that needs a stored entry the pattern lacks must fail rather than
    # decompress into the stale coloring.
    for detector in (TracerSparsityDetector(), TracerLocalSparsityDetector())
        backend = AutoSparse(
            AutoForwardDiff(); sparsity_detector = detector,
            coloring_algorithm = GreedyColoringAlgorithm()
        )
        alg = DFBDF(autodiff = backend, linsolve = LUFactorization())
        integrator = init(switching_prob, alg; solve_kwargs...)
        integrator.t = 0.75
        nlsolver = integrator.cache.nlsolver
        @test_throws DimensionMismatch OrdinaryDiffEqDifferentiation.calc_W!(
            nlsolver.cache.W, integrator, nlsolver, integrator.cache, 0.01, false, (true, true)
        )
        @test_throws DimensionMismatch solve(switching_prob, alg; abstol = 1.0e-8, reltol = 1.0e-8, solve_kwargs...)
    end

    # A known union is fixed by the caller, so it keeps the sparse path and stays correct.
    known = AutoSparse(
        AutoForwardDiff(); sparsity_detector = ADTypes.KnownJacobianSparsityDetector(union_pattern),
        coloring_algorithm = GreedyColoringAlgorithm()
    )
    alg = DFBDF(autodiff = known, linsolve = LUFactorization())
    sol = solve(switching_prob, alg; abstol = 1.0e-8, reltol = 1.0e-8, solve_kwargs...)
    @test successful_retcode(sol)
    @test sol.t[end] == 1.0
    @test maximum(abs, sol.u[end] - exact(1.0)) < 1.0e-6
    integrator = init(switching_prob, alg; abstol = 1.0e-8, reltol = 1.0e-8, solve_kwargs...)
    solve!(integrator)
    nlsolver = integrator.cache.nlsolver
    OrdinaryDiffEqDifferentiation.calc_W!(
        nlsolver.cache.W, integrator, nlsolver, integrator.cache, 0.01, false, (true, true)
    )
    @test Matrix(nlsolver.cache.J) ≈ J_u_late atol = 1.0e-8
    @test Matrix(nlsolver.cache.dae_jacobians.J_du) ≈ J_du atol = 1.0e-8
end

@testset "Explicit sparse jac_prototype with the default NoSparsityDetector" begin
    # The default detector detects nothing, so the prototype supplies the coloring.
    decoupled!(r, du, u, p, t) = (r .= du .+ u; nothing)
    backend = AutoSparse(AutoForwardDiff(); coloring_algorithm = GreedyColoringAlgorithm())
    for prototype in (sparse([1, 2], [1, 2], [1.0, 1.0]), sparse(ones(2, 2)))
        f = DAEFunction(decoupled!; jac_prototype = prototype)
        prob = DAEProblem(f, [-1.0, -1.0], [1.0, 1.0], (0.0, 0.5); differential_vars = [true, true])
        sol = solve(
            prob, DFBDF(autodiff = backend); initializealg = NoInit(),
            abstol = 1.0e-10, reltol = 1.0e-10
        )
        @test successful_retcode(sol)
        @test sol.u[end] ≈ fill(exp(-0.5), 2) rtol = 1.0e-6
    end
end

@testset "Tracer-detected sparsity on valid DAEs" begin
    local_sparse = AutoSparse(
        AutoForwardDiff(); sparsity_detector = TracerLocalSparsityDetector(),
        coloring_algorithm = GreedyColoringAlgorithm()
    )
    # The residual is only defined for du[1] > 1/2, so nothing may be traced at du = 0.
    function sqrt_residual!(r, du, u, p, t)
        r[1] = sqrt(du[1] - 0.5) - sqrt(0.5)
        r[2] = du[2] - 1
        return nothing
    end
    prototype = sparse([1, 2], [1, 2], [1.0, 1.0])
    prob = DAEProblem(
        DAEFunction(sqrt_residual!; jac_prototype = prototype), [1.0, 1.0], [1.0, 1.0],
        (0.0, 0.5); differential_vars = [true, true]
    )
    sol = solve(
        prob, DNordsieckBDF(autodiff = local_sparse); initializealg = NoInit(), dt = 0.125
    )
    @test successful_retcode(sol)
    @test sol.u[end] ≈ [1.5, 1.5] rtol = 1.0e-6

    # The construction trace samples du[2] = 8, where the branch is constant, so the
    # automatic storage misses (1, 2); at the real du[2] = -1 the partial is -2.
    function incomplete!(r, du, u, p, t)
        r[1] = du[1] + u[1] + (du[2] < 8 ? du[2]^2 - 1 : 0.0)
        r[2] = du[2] + u[2]
        return nothing
    end
    incomplete_prob = DAEProblem(
        incomplete!, [-1.0, -1.0], [1.0, 1.0], (0.0, 0.5); differential_vars = [true, true]
    )
    @test_throws "does not hold" solve(
        incomplete_prob, DFBDF(autodiff = local_sparse); initializealg = NoInit(), dt = 0.125
    )
end

# A detector may return a sparse Boolean pattern with stored `false` entries.
struct StoredFalseDetector <: ADTypes.AbstractSparsityDetector end
const stored_false_pattern = sparse([1, 2, 1], [1, 2, 2], [true, true, false], 2, 2)
ADTypes.jacobian_sparsity(f, x, ::StoredFalseDetector) = stored_false_pattern
ADTypes.jacobian_sparsity(f!, y, x, ::StoredFalseDetector) = stored_false_pattern

@testset "Stored false entries of a detected pattern detect nothing" begin
    decoupled!(r, du, u, p, t) = (r .= du .+ u; nothing)
    backend = AutoSparse(
        AutoForwardDiff(); sparsity_detector = StoredFalseDetector(),
        coloring_algorithm = GreedyColoringAlgorithm()
    )
    f = DAEFunction(decoupled!; jac_prototype = sparse([1, 2], [1, 2], [1.0, 1.0]))
    prob = DAEProblem(f, [-1.0, -1.0], [1.0, 1.0], (0.0, 0.5); differential_vars = [true, true])
    sol = solve(prob, DFBDF(autodiff = backend); initializealg = NoInit(), abstol = 1.0e-10, reltol = 1.0e-10)
    @test successful_retcode(sol)
    @test sol.u[end] ≈ fill(exp(-0.5), 2) rtol = 1.0e-6
end

@testset "A residual the tracer cannot trace at a refresh point still solves" begin
    # The global tracer fails on `u[1] < 0` once t >= 0.5, but the prototype covers
    # every derivative, so the solve must not fail because validation could not trace.
    function untraceable!(r, du, u, p, t)
        r[1] = du[1] + u[1] + (t >= 0.5 && u[1] < 0 ? u[1]^2 : 0.0)
        r[2] = du[2] + u[2]
        return nothing
    end
    backend = AutoSparse(
        AutoForwardDiff(); sparsity_detector = TracerSparsityDetector(),
        coloring_algorithm = GreedyColoringAlgorithm()
    )
    f = DAEFunction(untraceable!; jac_prototype = sparse([1, 2], [1, 2], [1.0, 1.0]))
    prob = DAEProblem(f, [-1.0, -1.0], [1.0, 1.0], (0.0, 1.0); differential_vars = [true, true])
    sol = solve(
        prob, DFBDF(autodiff = backend); initializealg = NoInit(),
        abstol = 1.0e-10, reltol = 1.0e-10, tstops = [0.5]
    )
    @test successful_retcode(sol)
    @test sol.u[end] ≈ fill(exp(-1.0), 2) rtol = 1.0e-6

    # Force a fresh Jacobian at a point the global tracer cannot trace.
    integrator = init(prob, DFBDF(autodiff = backend); initializealg = NoInit())
    integrator.t = 0.75
    integrator.uprev .= -1.0
    nlsolver = integrator.cache.nlsolver
    @test OrdinaryDiffEqDifferentiation.calc_W!(
        nlsolver.cache.W, integrator, nlsolver, integrator.cache, 0.01, false, (true, true)
    ) == (true, true)
    @test Matrix(nlsolver.cache.J) ≈ [-1.0 0.0; 0.0 1.0] atol = 1.0e-8

    # A pattern supplied as `sparsity` is user-owned storage too.
    f_sparsity = DAEFunction(untraceable!; sparsity = sparse([1, 2], [1, 2], [true, true]))
    prob_sparsity = DAEProblem(
        f_sparsity, [-1.0, -1.0], [1.0, 1.0], (0.0, 1.0); differential_vars = [true, true]
    )
    integrator = init(prob_sparsity, DFBDF(autodiff = backend); initializealg = NoInit())
    integrator.t = 0.75
    integrator.uprev .= -1.0
    nlsolver = integrator.cache.nlsolver
    @test OrdinaryDiffEqDifferentiation.calc_W!(
        nlsolver.cache.W, integrator, nlsolver, integrator.cache, 0.01, false, (true, true)
    ) == (true, true)
end

@testset "Untraceable refresh with automatically derived storage fails loudly" begin
    # No prototype: the initial trace (t = 0) stores the diagonal only. At t >= 0.5 the
    # tracer cannot evaluate `u[1] < 0`, and an unvalidated diagonal coloring would fold the
    # branch's dF1/du2 into the diagonal.
    function branching!(r, du, u, p, t)
        r[1] = du[1] + u[1] + (t >= 0.5 && u[1] < 0 ? u[2] : 0.0)
        r[2] = du[2] + u[2]
        return nothing
    end
    backend = AutoSparse(
        AutoForwardDiff(); sparsity_detector = TracerSparsityDetector(),
        coloring_algorithm = GreedyColoringAlgorithm()
    )
    prob = DAEProblem(
        branching!, [-1.0, -1.0], [1.0, 1.0], (0.0, 1.0); differential_vars = [true, true]
    )
    integrator = init(prob, DFBDF(autodiff = backend); initializealg = NoInit())
    integrator.t = 0.75
    integrator.uprev .= -1.0
    nlsolver = integrator.cache.nlsolver
    @test_throws "cannot trace the residual" OrdinaryDiffEqDifferentiation.calc_W!(
        nlsolver.cache.W, integrator, nlsolver, integrator.cache, 0.01, false, (true, true)
    )
end
