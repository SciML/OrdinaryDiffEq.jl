using DelayDiffEq, OrdinaryDiffEqBDF, Test

# With dt > lag every step overlaps the solution end, so DelayDiffEq recomputes
# each step at the same t in its fixed-point loop. u = c + t solves
# u'(t) = (u(t - τ) + τ) / (c + t) exactly, and BDF integrates linear solutions
# exactly. Starting away from t = 0 keeps the zero-filled history slots distinct
# from the real nodes, so reading an unfilled slot shows up as a wrong value.
@testset "FBDF time filter with fixed steps, t0 = $t0" for (t0, c) in
    ((0.0, 1.0), (-1.0, 2.0), (10.0, 1.0))
    τ = 0.125
    h(p, t) = [c + t]
    f(u, h, p, t) = (h(p, t - τ) .+ τ) ./ (c + t)
    f!(du, u, h, p, t) = (du .= f(u, h, p, t))
    h_scalar(p, t) = c + t
    f_scalar(u, h, p, t) = (h(p, t - τ) + τ) / (c + t)

    for (name, rhs, u0, history) in (
            ("in-place", f!, [c + t0], h),
            ("out-of-place array", f, [c + t0], h),
            ("out-of-place scalar", f_scalar, c + t0, h_scalar),
        )
        @testset "$name" begin
            prob = DDEProblem(rhs, u0, history, (t0, t0 + 1); constant_lags = [τ])
            sol = solve(
                prob, MethodOfSteps(FBDF(; time_filter = true)); dt = 0.2,
                adaptive = false
            )
            @test sol.retcode == ReturnCode.Success
            @test sol.t[end] == t0 + 1
            @test sol.stats.nfpiter > 0
            @test all(
                isapprox(only(u), c + t; atol = 1.0e-8, rtol = 1.0e-8)
                    for (u, t) in zip(sol.u, sol.t)
            )
        end
    end
end
