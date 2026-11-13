using DelayDiffEq, OrdinaryDiffEqSSPRK, StochasticDiffEqLowOrder
using SciMLBase: ReturnCode
using Test

# These caches only define `perform_step!(integrator, cache)`, so the fixed-point
# recompute must not pass `repeat_step`.

lag(u, p, t) = 1 + 0.1 * abs(first(u))
f(u, h, p, t) = -h(p, t - lag(u, p, t))
f!(du, u, h, p, t) = (du[1] = -h(p, t - lag(u, p, t))[1]; nothing)
g(u, h, p, t) = 0.1
g!(du, u, h, p, t) = (du[1] = 0.1; nothing)

@testset "ODE: KYKSSPRK42" begin
    prob = DDEProblem(f, 1.0, (p, t) -> 1.0, (0.0, 5.0); dependent_lags = (lag,))
    sol = solve(prob, MethodOfSteps(KYKSSPRK42()); dt = 1.5)
    @test sol.retcode == ReturnCode.Success
    @test sol.stats.nfpiter > 0
end

@testset "SDE: LambaEM, out-of-place" begin
    prob = SDDEProblem(f, g, 1.0, (p, t) -> 1.0, (0.0, 5.0); dependent_lags = (lag,))
    sol = solve(prob, MethodOfSteps(LambaEM()); seed = 1)
    @test sol.retcode == ReturnCode.Success
    @test sol.stats.nfpiter > 0
end

@testset "SDE: LambaEM, in-place" begin
    prob = SDDEProblem(f!, g!, [1.0], (p, t) -> [1.0], (0.0, 5.0); dependent_lags = (lag,))
    sol = solve(prob, MethodOfSteps(LambaEM()); seed = 1)
    @test sol.retcode == ReturnCode.Success
    @test sol.stats.nfpiter > 0
end
