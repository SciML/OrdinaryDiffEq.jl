using DelayDiffEq, OrdinaryDiffEqBDF, Test

@testset "FBDF overlapping adaptive steps" begin
    h(p, t) = [1 + t]
    f(u, h, p, t) = (h(p, t - 0.125) .+ 0.125) ./ (1 + t)
    f!(du, u, h, p, t) = (du .= f(u, h, p, t))
    h_scalar(p, t) = 1 + t
    f_scalar(u, h, p, t) = (h(p, t - 0.125) + 0.125) / (1 + t)

    for (name, rhs, u0, history) in (
            ("in-place", f!, [1.0], h),
            ("out-of-place array", f, [1.0], h),
            ("out-of-place scalar", f_scalar, 1.0, h_scalar),
        )
        @testset "$name" begin
            prob = DDEProblem(rhs, u0, history, (0.0, 1.0); constant_lags = [0.125])
            sol = solve(
                prob, MethodOfSteps(FBDF()); dt = 0.2, adaptive = true,
                abstol = 1.0e-8, reltol = 1.0e-8
            )
            @test sol.retcode == ReturnCode.Success
            @test sol.t[end] == 1.0
            @test sol.stats.nfpiter > 0
            @test all(
                isapprox(only(u), 1 + t; atol = 1.0e-8, rtol = 1.0e-8)
                    for (u, t) in zip(sol.u, sol.t)
            )
        end
    end
end
