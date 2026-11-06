using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, Test
using SciMLBase: ReturnCode

@testset "Float32 short tspan: no oversized tstop snap" begin
    f(u, p, t) = p .* u
    t0 = 1.5f0
    tf = 1.50001f0
    p = Float32[-6.0f6]
    s = solve(ODEProblem(f, Float32[1], (t0, tf), p), Tsit5(); dt = 8 * eps(t0))
    @test s.retcode == ReturnCode.Success
    @test s.stats.naccept > 1
    @test s.u[end][1] < 1.0f-5
end

@testset "Snap reaches tstop: u' = 1 covers full span" begin
    f(u, p, t) = one(u)
    t0 = 1.5f0
    span = 8 * eps(t0)
    tf = t0 + span
    dt = 6 * eps(t0)
    s = solve(ODEProblem(f, zero(t0), (t0, tf)), Tsit5(); dt = dt, adaptive = false)
    @test s.retcode == ReturnCode.Success
    @test s.t[end] == tf
    @test s.u[end] ≈ (tf - t0)

    s_rev = solve(
        ODEProblem(f, zero(t0), (tf, t0)), Tsit5();
        dt = -dt, adaptive = false
    )
    @test s_rev.retcode == ReturnCode.Success
    @test s_rev.t[end] == t0
    @test s_rev.u[end] ≈ (t0 - tf)

    t0_64 = 1.0e10
    span_64 = 150 * eps(t0_64)
    tf_64 = t0_64 + span_64
    dt_64 = 110 * eps(t0_64)
    s64 = solve(
        ODEProblem(f, 0.0, (t0_64, tf_64)), Tsit5();
        dt = dt_64, adaptive = false
    )
    @test s64.retcode == ReturnCode.Success
    @test s64.t[end] == tf_64
    @test s64.u[end] ≈ (tf_64 - t0_64)
end

@testset "Adaptive snap reaches tstop: u' = 1 covers full span" begin
    f(u, p, t) = one(u)
    t0 = 1.5f0
    span = 8 * eps(t0)
    tf = t0 + span
    dt = 6 * eps(t0)
    s = solve(ODEProblem(f, zero(t0), (t0, tf)), Tsit5(); dt = dt)
    @test s.retcode == ReturnCode.Success
    @test s.u[end] ≈ (tf - t0)

    t0_64 = 1.0e10
    span_64 = 150 * eps(t0_64)
    tf_64 = t0_64 + span_64
    dt_64 = 110 * eps(t0_64)
    s64 = solve(ODEProblem(f, 0.0, (t0_64, tf_64)), Tsit5(); dt = dt_64)
    @test s64.retcode == ReturnCode.Success
    @test s64.u[end] ≈ (tf_64 - t0_64)
end

@testset "Float32/Float64 tstops are hit without collapsing steps" begin
    f(u, p, t) = -u
    counts = Int[]
    for (T, tspan, tstops) in (
            (Float64, (0.0, 1.0), [0.25, 0.5, 0.75]),
            (Float32, (0.0f0, 1.0f0), Float32[0.25, 0.5, 0.75]),
        )
        prob = ODEProblem(f, T(1), tspan)
        s = solve(prob, Tsit5(); tstops = tstops, save_everystep = true)
        @test s.retcode == ReturnCode.Success
        @test all(t -> t in s.t, tstops)
        @test s.t[end] == tspan[2]
        push!(counts, s.stats.naccept)
    end
    @test counts[1] == counts[2]
    @test counts[1] >= length([0.25, 0.5, 0.75])
end

@testset "First-step rejection near a tstop" begin
    t0 = 1.5f0
    span = 8 * eps(t0)
    for factor in (1.2f0, 1.4f0, 1.6f0), direction in (1, -1)
        tspan = direction > 0 ? (t0, t0 + span) : (t0 + span, t0)
        prob = ODEProblem((u, p, t) -> -p * u, 1.0f0, tspan, direction * factor / span)
        sol = solve(prob, Tsit5(); dt = direction * 6 * eps(t0), maxiters = 30)
        @test sol.retcode == ReturnCode.Success
        @test sol.stats.nreject > 0
        @test sol.t[end] == tspan[2]
        @test sol.u[end] ≈ exp(-factor) rtol = 2.0f-3
    end
end

@testset "Snapped steps respect dtmax" begin
    t0 = 1.5f0
    tf = t0 + 8 * eps(t0)
    dtmax = 6 * eps(t0)
    for adaptive in (false, true), direction in (1, -1), dt in (dtmax, 10 * eps(t0))
        tspan = direction > 0 ? (t0, tf) : (tf, t0)
        prob = ODEProblem((u, p, t) -> one(u), 0.0f0, tspan)
        sol = solve(prob, Tsit5(); adaptive, dt = direction * dt, dtmax)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] == tspan[2]
        @test maximum(abs, diff(sol.t)) <= dtmax
        @test sol.u[end] ≈ tspan[2] - tspan[1]
    end
end
