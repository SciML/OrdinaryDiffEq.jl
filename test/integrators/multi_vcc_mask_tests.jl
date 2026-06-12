using OrdinaryDiffEq, Test

# Regression tests for https://github.com/SciML/OrdinaryDiffEq.jl/issues/3673
# When multiple VectorContinuousCallbacks share a CallbackSet, the `events` mask
# passed to each callback's `affect!` must:
#   - have length equal to that callback's `len` (not the max across all VCCs)
#   - reflect that callback's own triggered events (not those of a sibling VCC
#     processed later in `find_first_continuous_callback`)

@testset "Multiple VectorContinuousCallbacks: mask correctness (#3673)" begin
    prob = ODEProblem((du, u, p, t) -> (du[1] = 1.0), [0.0], (0.0, 10.0))

    # VCC #1: two conditions, upcrossings at u = 3 and u = 5
    function vcc1_cond(out, u, _, _)
        out[1] = u[1] - 3.0
        out[2] = u[1] - 5.0
    end
    fired1 = Tuple{Float64, Vector{Int8}}[]
    vcc1_affect!(integ, events) = push!(fired1, (integ.t, collect(events)))
    cb1 = VectorContinuousCallback(vcc1_cond, vcc1_affect!, 2)

    # VCC #2: one condition, upcrossing at u = 8
    function vcc2_cond(out, u, _, _)
        out[1] = u[1] - 8.0
    end
    fired2 = Tuple{Float64, Vector{Int8}}[]
    vcc2_affect!(integ, events) = push!(fired2, (integ.t, collect(events)))
    cb2 = VectorContinuousCallback(vcc2_cond, vcc2_affect!, 1)

    # VCC #3: three conditions, upcrossings at u = 4, u = 6, u = 8
    function vcc3_cond(out, u, _, _)
        out[1] = u[1] - 4.0
        out[2] = u[1] - 6.0
        out[3] = u[1] - 8.0
    end
    fired3 = Tuple{Float64, Vector{Int8}}[]
    vcc3_affect!(integ, events) = push!(fired3, (integ.t, collect(events)))
    cb3 = VectorContinuousCallback(vcc3_cond, vcc3_affect!, 3)

    @testset "cb1 (len=2) + cb2 (len=1)" begin
        empty!(fired1); empty!(fired2)
        solve(prob, Tsit5(), callback = CallbackSet(cb1, cb2))

        @test length(fired1) == 2
        @test length(fired2) == 1

        # All masks must match each callback's own `len`
        @test all(length(m) == 2 for (_, m) in fired1)
        @test all(length(m) == 1 for (_, m) in fired2)

        # cb1's first firing is at u=3 (idx 1), second at u=5 (idx 2)
        @test fired1[1][1] ≈ 3.0 atol=1e-10
        @test fired1[1][2] == Int8[1, 0]
        @test fired1[2][1] ≈ 5.0 atol=1e-10
        @test fired1[2][2] == Int8[0, 1]

        # cb2 fires at u=8 (idx 1)
        @test fired2[1][1] ≈ 8.0 atol=1e-10
        @test fired2[1][2] == Int8[1]
    end

    @testset "cb1 (len=2) + cb3 (len=3)" begin
        empty!(fired1); empty!(fired3)
        solve(prob, Tsit5(), callback = CallbackSet(cb1, cb3))

        @test length(fired1) == 2
        @test length(fired3) >= 3

        @test all(length(m) == 2 for (_, m) in fired1)
        @test all(length(m) == 3 for (_, m) in fired3)

        @test fired1[1][1] ≈ 3.0 atol=1e-10
        @test fired1[1][2] == Int8[1, 0]
        @test fired1[2][1] ≈ 5.0 atol=1e-10
        @test fired1[2][2] == Int8[0, 1]

        # cb3's idx 1 fires at u=4, idx 2 at u=6, idx 3 at u=8
        idx4 = findfirst(((t, _),) -> isapprox(t, 4.0; atol=1e-10), fired3)
        idx6 = findfirst(((t, _),) -> isapprox(t, 6.0; atol=1e-10), fired3)
        idx8 = findfirst(((t, _),) -> isapprox(t, 8.0; atol=1e-10), fired3)
        @test idx4 !== nothing
        @test idx6 !== nothing
        @test idx8 !== nothing
        @test fired3[idx4][2] == Int8[1, 0, 0]
        @test fired3[idx6][2] == Int8[0, 1, 0]
        @test fired3[idx8][2] == Int8[0, 0, 1]
    end

    # Reverse the order so the larger VCC is processed first; the shared
    # `callback_cache.simultaneous_events` buffer is sized to max(len), so the
    # bug previously surfaced regardless of order. Re-verify masks.
    @testset "cb3 (len=3) + cb1 (len=2) — reverse order" begin
        empty!(fired1); empty!(fired3)
        solve(prob, Tsit5(), callback = CallbackSet(cb3, cb1))

        @test all(length(m) == 2 for (_, m) in fired1)
        @test all(length(m) == 3 for (_, m) in fired3)
        @test fired1[1][2] == Int8[1, 0]
        @test fired1[2][2] == Int8[0, 1]
    end
end
