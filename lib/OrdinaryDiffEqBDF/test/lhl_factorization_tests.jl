using OrdinaryDiffEqBDF, LinearSolve, LinearAlgebra, SparseArrays, Test
using SciMLOperators: WOperator, jacobian_stale

# `LHLFactorization` needs W left split as J plus a scalar shift so that a new dtgamma
# costs O(n²).  These tests pin down that the split form is what the integrator builds,
# that it produces the same trajectory as an assembled W, and that unsupported mass
# matrices are rejected rather than silently densified.

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

@testset "sparse Jacobian is kept split" begin
    A = sparse([1, 2, 3, 2, 3, 1], [1, 2, 3, 1, 2, 3], [-4.0, -5.0, -6.0, 1.0, 1.0, 1.0], 3, 3)
    f! = (du, u, p, t) -> mul!(du, A, u)
    jac! = (J, u, p, t) -> copyto!(J, A)
    prob = ODEProblem(ODEFunction(f!; jac = jac!, jac_prototype = copy(A)), ones(3), (0.0, 0.1))
    integ = init(prob, FBDF(linsolve = LHLFactorization()); abstol = 1.0e-10, reltol = 1.0e-10)

    @test integ.cache.nlsolver.cache.W isa WOperator
    @test integ.cache.nlsolver.cache.W.J isa SparseMatrixCSC
    sol = solve!(integ)
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end] ≈ exp(0.1 * Matrix(A)) * ones(3) rtol = 1.0e-7
end

@testset "unsupported mass matrix is rejected" begin
    fm = ODEFunction((du, u, p, t) -> (du .= -u); mass_matrix = Diagonal([1.0, 2.0, 3.0]))
    mprob = ODEProblem(fm, ones(3), (0.0, 1.0))
    @test_throws ArgumentError solve(mprob, NordsieckBDF(linsolve = LHLFactorization()))
end
