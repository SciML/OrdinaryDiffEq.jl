using Test, DelayDiffEq, OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqHighOrderRK,
    OrdinaryDiffEqRosenbrock

# Regression test for SciML/OrdinaryDiffEq.jl#3143.
#
# Several methods, e.g. DP5, DP8 and the Rosenbrock methods, only compute the
# coefficients of their dense output when `integrator.opts.calck` is true. `calck` used to default to `dense`,
# which is false with `saveat` or `save_everystep = false`, so the history function
# interpolated stale coefficients and the solution lost several orders of accuracy.
# The history needs the interpolation data regardless of what is saved: saving less
# must not change the steps or the solution values.

f_3143!(du, u, h, p, t) = (du[1] = u[1] * (1 - h(p, t - p)[1]); nothing)
f_3143(u, h, p, t) = [u[1] * (1 - h(p, t - p)[1])]
h_3143(p, t) = [1.2]

@testset "issue #3143: history without dense output ($(nameof(typeof(alg))), iip = $iip)" for alg in (
            DP5(), DP8(), Rodas4(),
        ),
        iip in (true, false)

    prob = DDEProblem(
        iip ? f_3143! : f_3143, [1.2], h_3143, (0.0, 20.0), 1.0; constant_lags = [1.0]
    )
    ddealg = MethodOfSteps(alg)
    ts = 0.0:0.5:20.0
    sol_dense = solve(prob, ddealg)

    sol_saveat = solve(prob, ddealg; saveat = ts)
    @test sol_saveat.stats.naccept == sol_dense.stats.naccept
    @test sol_saveat.t == collect(ts)
    @test all(isapprox(sol_saveat.u[i], sol_dense(t); atol = 1.0e-10) for (i, t) in enumerate(ts))

    sol_last = solve(prob, ddealg; save_everystep = false)
    @test sol_last.stats.naccept == sol_dense.stats.naccept
    @test isapprox(sol_last.u[end], sol_dense.u[end]; atol = 1.0e-10)
end
