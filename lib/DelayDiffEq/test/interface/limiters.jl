using DelayDiffEq
using OrdinaryDiffEqSSPRK
using Test

const LIMITER_CALLS = Ref(0)
limiter!(u, integrator, p, t) = (LIMITER_CALLS[] += 1)

f(du, u, h, p, t) = (du[1] = -0.5 * u[1])
hist(p, t) = [1.0]
prob = DDEProblem(f, [1.0], hist, (0.0, 3.0))

@testset "solve-level step_limiter keyword" begin
    LIMITER_CALLS[] = 0
    sol = solve(prob, MethodOfSteps(SSPRK43()), dt = 0.1; step_limiter = limiter!)
    @test LIMITER_CALLS[] > 0
    @test sol.stats.naccept + sol.stats.nreject == LIMITER_CALLS[]
end

@testset "deprecated per-algorithm step_limiter! field is honored" begin
    LIMITER_CALLS[] = 0
    sol = solve(prob, MethodOfSteps(SSPRK43(; step_limiter! = limiter!)), dt = 0.1)
    @test LIMITER_CALLS[] > 0
    @test sol.stats.naccept + sol.stats.nreject == LIMITER_CALLS[]
end
