using OrdinaryDiffEqBDF, OrdinaryDiffEqFIRK, LinearSolve, LinearAlgebra, SparseArrays, Test
using SciMLOperators: WOperator, jacobian_stale

# `LHLFactorization` needs W left split as J plus a scalar shift so that a new dtgamma
# costs O(n²).  These tests pin down that the split form is what the integrator builds,
# that it produces the same trajectory as an assembled W, and that the combinations the
# reduction cannot serve are rejected rather than silently densified.

function rober_like(n)
    # A stiff, strongly coupled system with a dense Jacobian and no sparsity to exploit.
    f = function (du, u, p, t)
        @inbounds for k in 0:(n - 1)
            i = 3k
            du[i + 1] = -0.04u[i + 1] + 1.0e4 * u[i + 2] * u[i + 3]
            du[i + 2] = 0.04u[i + 1] - 1.0e4 * u[i + 2] * u[i + 3] - 3.0e7 * u[i + 2]^2
            du[i + 3] = 3.0e7 * u[i + 2]^2
        end
        return nothing
    end
    u0 = repeat([1.0, 0.0, 0.0], n)
    return ODEProblem(ODEFunction(f), u0, (0.0, 1.0e4))
end

const PROB = rober_like(4)

@testset "W is kept split" begin
    integ = init(PROB, NordsieckBDF(linsolve = LHLFactorization()); abstol = 1.0e-8, reltol = 1.0e-8)
    W = integ.cache.nlsolver.cache.W
    @test W isa WOperator
    @test W.J === integ.cache.nlsolver.cache.J
    integ2 = init(PROB, NordsieckBDF(linsolve = LUFactorization()); abstol = 1.0e-8, reltol = 1.0e-8)
    @test integ2.cache.nlsolver.cache.W isa Matrix
end

@testset "same trajectory as an assembled W" begin
    for alg in (NordsieckBDF, FBDF, QNDF)
        lu = solve(PROB, alg(linsolve = LUFactorization()); abstol = 1.0e-10, reltol = 1.0e-10)
        lhl = solve(PROB, alg(linsolve = LHLFactorization()); abstol = 1.0e-10, reltol = 1.0e-10)
        @test SciMLBase.successful_retcode(lhl)
        @test lhl.u[end] ≈ lu.u[end] rtol = 1.0e-6
    end
end

@testset "the reduction is reused across dtgamma changes" begin
    integ = init(PROB, NordsieckBDF(linsolve = LHLFactorization()); abstol = 1.0e-8, reltol = 1.0e-8)
    solve!(integ)
    # More W rebuilds than Jacobians is the whole point; each extra one is O(n²).
    @test integ.stats.nw > integ.stats.njacs
    lc = integ.cache.nlsolver.cache.linsolve.cacheval
    @test lc isa LinearSolve.LHLCache
    # The linear solve reduced and claimed the operator, so nothing is left outstanding,
    # and the reduction is stamped with the Jacobian it was actually taken of.
    @test lc.ws.reduced
    @test lc.jac === integ.cache.nlsolver.cache.J
    @test !jacobian_stale(integ.cache.nlsolver.cache.W)
end

@testset "Radau keeps both stage matrices split" begin
    # RadauIIA forms a real stage matrix and a complex-conjugate pair from one Jacobian —
    # exactly the workload the reduction exists for. Both stages stay split: the real one
    # as the WOperator `build_J_W` made, the complex pair as a WOperator whose `gamma`
    # slot is complex while its Jacobian — the *same object* — stays real, so inside the
    # linear solver one real reduction serves a real and a complex shift.
    integ = init(PROB, RadauIIA5(linsolve = LHLFactorization()); abstol = 1.0e-8, reltol = 1.0e-8)
    (; J, W1, W2) = integ.cache
    @test W1 isa WOperator
    @test W2 isa WOperator
    @test W1.J === J
    @test W2.J === J
    @test typeof(W2.gamma) <: Complex

    for i in 1:6
        step!(integ)
    end
    ws1 = integ.cache.linsolve1.cacheval.ws
    ws2 = integ.cache.linsolve2.cacheval.ws
    @test eltype(ws1.factors) === Float64   # both reductions stay real...
    @test eltype(ws2.factors) === Float64
    @test typeof(ws2.σ) <: Complex          # ...only the pair's shift is complex

    lu = solve(PROB, RadauIIA5(linsolve = LUFactorization()); abstol = 1.0e-10, reltol = 1.0e-10)
    for alg in (RadauIIA3, RadauIIA5, RadauIIA9, AdaptiveRadau)
        lhl = solve(PROB, alg(linsolve = LHLFactorization()); abstol = 1.0e-10, reltol = 1.0e-10)
        @test SciMLBase.successful_retcode(lhl)
        @test lhl.u[end] ≈ lu.u[end] rtol = 1.0e-6
    end

    # AdaptiveRadau holds its complex pairs in a *vector* of split Ws, one per conjugate
    # pair up to the maximum stage count, all viewing the one Jacobian; the in-place
    # refresh must reach every element of that vector, not just the scalar fields.
    ainteg = init(
        PROB, AdaptiveRadau(linsolve = LHLFactorization()); abstol = 1.0e-8, reltol = 1.0e-8
    )
    @test ainteg.cache.W1 isa WOperator
    @test ainteg.cache.W2 isa AbstractVector
    @test all(W -> W isa WOperator && W.J === ainteg.cache.J, ainteg.cache.W2)
    @test all(W -> typeof(W.gamma) <: Complex, ainteg.cache.W2)
end

@testset "Radau re-reduces when the Jacobian moves" begin
    # The stage matrices share one Jacobian that `firk_new_J!` rewrites in place; every
    # split W must be marked stale or the pair's cached reduction would keep answering
    # for the previous Jacobian. rober_like refreshes J many times over this span, so a
    # missed mark shows up as a wrong trajectory, not a subtle residual.
    integ = init(PROB, RadauIIA5(linsolve = LHLFactorization()); abstol = 1.0e-8, reltol = 1.0e-8)
    solve!(integ)
    @test SciMLBase.successful_retcode(integ.sol)
    @test integ.sol.stats.njacs > 1
end

@testset "unsupported combinations are rejected" begin
    f = ODEFunction((du, u, p, t) -> (du .= -u); jac_prototype = spzeros(3, 3))
    sprob = ODEProblem(f, ones(3), (0.0, 1.0))
    @test_throws ArgumentError solve(sprob, NordsieckBDF(linsolve = LHLFactorization()))

    fm = ODEFunction((du, u, p, t) -> (du .= -u); mass_matrix = Diagonal([1.0, 2.0, 3.0]))
    mprob = ODEProblem(fm, ones(3), (0.0, 1.0))
    @test_throws ArgumentError solve(mprob, NordsieckBDF(linsolve = LHLFactorization()))
end
