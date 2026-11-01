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
