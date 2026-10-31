using OrdinaryDiffEqSDIRK, OrdinaryDiffEqBDF, OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg, get_linear_cache
using NonlinearSolve: NewtonRaphson
using LinearSolve, LinearAlgebra, ADTypes, SciMLBase
using SciMLOperators: WOperator, jacobian_stale
using Test

# `LHLFactorization` reduces `J` once and absorbs each new `γΔt` in O(n²), which only works
# while `W` stays split as `J` plus a scalar shift. NonlinearSolveAlg used to hand the inner
# solver an assembled copy of that `W`, so the inner factorization saw a bare matrix and
# redid the whole reduction on every refresh.

function rober_like(n)
    f = function (du, u, p, t)
        @inbounds for k in 0:(n - 1)
            i = 3k
            du[i + 1] = -0.04u[i + 1] + 1.0e4 * u[i + 2] * u[i + 3]
            du[i + 2] = 0.04u[i + 1] - 1.0e4 * u[i + 2] * u[i + 3] - 3.0e7 * u[i + 2]^2
            du[i + 3] = 3.0e7 * u[i + 2]^2
        end
        return nothing
    end
    return ODEProblem(ODEFunction(f), repeat([1.0, 0.0, 0.0], n), (0.0, 1.0e4))
end

const PROB = rober_like(4)
nsa(; kw...) = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoForwardDiff()); kw...)

@testset "the inner solver receives the split W itself" begin
    integ = init(
        PROB, TRBDF2(linsolve = LHLFactorization(), nlsolve = nsa());
        abstol = 1.0e-8, reltol = 1.0e-8
    )
    c = integ.cache.nlsolver.cache
    @test c.W isa WOperator
    @test c.W.J === c.J
    lc = get_linear_cache(c.cache)
    @test lc.A === c.W
    @test c.linsolve === lc
    step!(integ)
    @test lc.cacheval isa LinearSolve.LHLCache
    @test lc.cacheval.jac === c.J

    integ_lu = init(
        PROB, TRBDF2(linsolve = LUFactorization(), nlsolve = nsa());
        abstol = 1.0e-8, reltol = 1.0e-8
    )
    @test integ_lu.cache.nlsolver.cache.W isa Matrix
end

@testset "same trajectory as an assembled W" begin
    for alg in (TRBDF2, KenCarp4, FBDF)
        lu = solve(
            PROB, alg(linsolve = LUFactorization(), nlsolve = nsa());
            abstol = 1.0e-10, reltol = 1.0e-10
        )
        lhl = solve(
            PROB, alg(linsolve = LHLFactorization(), nlsolve = nsa());
            abstol = 1.0e-10, reltol = 1.0e-10
        )
        @test SciMLBase.successful_retcode(lhl)
        @test lhl.u[end] ≈ lu.u[end] rtol = 1.0e-6
    end
end

@testset "the reduction is reused across γΔt changes" begin
    integ = init(
        PROB, TRBDF2(linsolve = LHLFactorization(), nlsolve = nsa());
        abstol = 1.0e-8, reltol = 1.0e-8
    )
    solve!(integ)
    @test integ.stats.nw > integ.stats.njacs
    lc = integ.cache.nlsolver.cache.linsolve.cacheval
    @test lc.ws.reduced
    @test lc.jac === integ.cache.nlsolver.cache.J
    @test !jacobian_stale(integ.cache.nlsolver.cache.W)
end

@testset "split W takes every γΔt change" begin
    A = [-20.0 1.0; 1.0 -20.0]
    f!(du, u, p, t) = mul!(du, A, u)
    prob = ODEProblem(f!, [1.0, 0.5], (0.0, 1.0))
    lhl(cutoff) = SDIRK2(
        linsolve = LHLFactorization(), nlsolve = nsa(new_W_dt_cutoff = cutoff)
    )
    loose = solve(prob, lhl(0.2); abstol = 1.0e-6, reltol = 1.0e-6)
    exact = solve(prob, lhl(0.0); abstol = 1.0e-6, reltol = 1.0e-6)
    @test loose.u[end] == exact.u[end]
    @test loose.stats.naccept == exact.stats.naccept
    @test loose.stats.nw == exact.stats.nw
end
