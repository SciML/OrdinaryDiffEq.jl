using OrdinaryDiffEqSDC
using SciMLBase
using Test

const L1 = -2.0
const L2 = 0.5

imex_oop() = SplitODEProblem(
    SplitFunction(
        (u, p, t) -> L1 * u, (u, p, t) -> L2 * u;
        analytic = (u0, p, t) -> u0 * exp((L1 + L2) * t)
    ),
    1.0, (0.0, 1.0)
)
imex_iip() = SplitODEProblem(
    SplitFunction(
        (du, u, p, t) -> (du .= L1 .* u), (du, u, p, t) -> (du .= L2 .* u);
        analytic = (u0, p, t) -> u0 .* exp((L1 + L2) * t)
    ),
    [1.0, 0.5], (0.0, 1.0)
)

function observed_order(prob, alg)
    nsteps = [8, 16, 32, 64]
    exact(t) = prob.f.analytic(prob.u0, prob.p, t)
    errs = map(nsteps) do n
        sol = solve(prob, alg; dt = 1 / n, adaptive = false)
        maximum(abs, sol.u[end] .- exact(1.0))
    end
    return log2(errs[end - 1] / errs[end])
end

@testset "IMEX order gains one per sweep" begin
    for (make, tag) in ((imex_oop, "oop"), (imex_iip, "iip")), K in 1:4
        alg = SDC(num_nodes = 3, num_sweeps = K, sweeper = SDCSweeper.BE)
        @testset "$tag K=$K" begin
            @test observed_order(make(), alg) ≈ K + 1 atol = 0.3
        end
    end
end

@testset "the explicit part is never differentiated" begin
    function f2_no_duals(u, p, t)
        eltype(u) <: AbstractFloat || error("f2 reached automatic differentiation")
        return L2 * u
    end
    function f2_no_duals!(du, u, p, t)
        eltype(u) <: AbstractFloat || error("f2 reached automatic differentiation")
        du .= L2 .* u
        return nothing
    end
    prob_oop = SplitODEProblem((u, p, t) -> L1 * u, f2_no_duals, 1.0, (0.0, 1.0))
    prob_iip = SplitODEProblem(
        (du, u, p, t) -> (du .= L1 .* u), f2_no_duals!, [1.0, 0.5], (0.0, 1.0)
    )
    for prob in (prob_oop, prob_iip)
        sol = solve(prob, SDC(num_sweeps = 3); dt = 0.1, adaptive = false)
        @test SciMLBase.successful_retcode(sol)
        @test sol.stats.nf2 > 0
    end
end

@testset "threaded IMEX needs a Picard explicit sweeper" begin
    prob = imex_iip()
    alg = SDC(sweeper = SDCSweeper.MIN_SR_S, threading = true)
    @test_throws ArgumentError init(prob, alg; dt = 0.1, adaptive = false)
    alg = SDC(
        sweeper = SDCSweeper.MIN_SR_S, explicit_sweeper = SDCSweeper.Picard,
        threading = true
    )
    @test SciMLBase.successful_retcode(solve(prob, alg; dt = 0.1, adaptive = false))
end
