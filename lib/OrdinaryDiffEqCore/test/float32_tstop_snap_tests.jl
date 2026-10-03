using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, Test
using SciMLBase: ReturnCode

@testset "Float32 short tspan: no oversized tstop snap" begin
    # Adaptive Tsit5 with an explicit small dt on a Float32 span shorter than
    # 100*eps(t) must take many steps, not teleport to tf after one accept.
    f(u, p, t) = p .* u
    t0 = 1.5f0
    tf = 1.50001f0
    p = Float32[-6.0f6]
    s = solve(ODEProblem(f, Float32[1], (t0, tf), p), Tsit5(); dt = 8 * eps(t0))
    @test s.retcode == ReturnCode.Success
    @test s.stats.naccept > 1
    # Exact decay is ~8e-27; Float32 roundoff leaves a small residual, not O(0.1).
    @test s.u[end][1] < 1.0f-5
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
    # Same problem in Float32/Float64 should take the same number of accepts.
    @test counts[1] == counts[2]
    @test counts[1] >= length([0.25, 0.5, 0.75])
end
