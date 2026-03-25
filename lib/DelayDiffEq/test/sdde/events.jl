using DelayDiffEq
using StochasticDiffEqLowOrder
using Test

@testset "SDDE discrete callback" begin
    f(u, h, p, t) = 0.3 * h(p, t - p)
    g(u, h, p, t) = 0.1
    h(p, t) = 0.0
    prob = SDDEProblem(f, g, h, (0.0, 4.0), 0.5; constant_lags = (0.5,))

    # event at `t = 3`
    cb = DiscreteCallback(
        (u, t, integrator) -> t == 3,
        integrator -> (integrator.u = -integrator.u)
    )

    sol = solve(prob, MethodOfSteps(RKMil()), dt = 0.01, adaptive = false,
        tstops = (3,), callback = cb)
    ts = findall(x -> x == 3, sol.t)
    @test length(ts) == 2
    @test sol.u[ts[1]] == -sol.u[ts[2]]
end
