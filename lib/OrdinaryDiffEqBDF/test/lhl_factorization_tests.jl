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

@testset "FIRK refuses a reduction it cannot reuse" begin
    # RadauIIA assembles a real and a complex stage matrix per step, so an LHL reduction
    # would be re-taken for each one rather than reused: measured 5.2x slower than LU on a
    # 128-unknown Brusselator, at identical step counts and error, with 279 stage-matrix
    # updates against 3 Jacobians. Refuse rather than silently pessimize; reusing one
    # reduction across both is SciML/OrdinaryDiffEq.jl#4281.
    @test_throws ArgumentError solve(PROB, RadauIIA5(linsolve = LHLFactorization()))
    # The supported solvers are unaffected.
    lu = solve(PROB, RadauIIA5(linsolve = LUFactorization()); abstol = 1.0e-10, reltol = 1.0e-10)
    @test SciMLBase.successful_retcode(lu)
end

@testset "unsupported combinations are rejected" begin
    f = ODEFunction((du, u, p, t) -> (du .= -u); jac_prototype = spzeros(3, 3))
    sprob = ODEProblem(f, ones(3), (0.0, 1.0))
    @test_throws ArgumentError solve(sprob, NordsieckBDF(linsolve = LHLFactorization()))

    fm = ODEFunction((du, u, p, t) -> (du .= -u); mass_matrix = Diagonal([1.0, 2.0, 3.0]))
    mprob = ODEProblem(fm, ones(3), (0.0, 1.0))
    @test_throws ArgumentError solve(mprob, NordsieckBDF(linsolve = LHLFactorization()))
end
