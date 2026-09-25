# The scalar `ode_interpolation`/`ode_interpolation!` paths warm-start their
# interval search from a `TsSearchHint` stored in `InterpolationData`. These
# tests pin two properties:
#   1. Correctness: scalar interpolation matches the (hint-free) vector-tvals
#      path exactly, for ascending / descending / random access orders, on
#      both forward-time and reverse-time solutions, for both continuities.
#   2. The warm start actually engages: after a scalar evaluation the hint
#      holds the last hit instead of its initial value.
using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, SciMLBase, Random, Test

# Damped rotation: well-conditioned in both time directions.
function rot!(du, u, p, t)
    du[1] = -0.1 * u[1] - u[2]
    du[2] = u[1] - 0.1 * u[2]
    return nothing
end

@testset "scalar interpolation hint" begin
    tq = collect(range(0.01, 9.99, length = 499))
    orders = [
        ("ascending", tq),
        ("descending", reverse(tq)),
        ("random", tq[randperm(Xoshiro(1), length(tq))]),
    ]
    for tspan in ((0.0, 10.0), (10.0, 0.0))
        sol = solve(
            ODEProblem(rot!, [1.0, 0.0], tspan), Tsit5();
            abstol = 1.0e-10, reltol = 1.0e-10, dense = true
        )
        # reference through the (unchanged) vector-tvals path, per continuity
        ref = Dict(
            cont => Dict(zip(tq, sol(tq; continuity = cont).u))
                for cont in (:left, :right)
        )
        @testset "tspan=$tspan, $name, continuity=$cont" for (name, tvals) in orders,
                cont in (:left, :right)

            @test all(sol(t; continuity = cont) == ref[cont][t] for t in tvals)
            # in-place form
            y = zeros(2)
            ok = all(tvals) do t
                sol(y, t; continuity = cont)
                y == ref[cont][t]
            end
            @test ok
        end
        # the hint engages: a scalar evaluation moves it off the initial index
        hint = sol.interp.ts_hint
        hint.idx_prev = 1
        sol(tq[end ÷ 2])
        @test hint.idx_prev != 1
        # a fresh hint has no previous hit and guesses `lastindex(ts)` at query
        # time: the grid grows after construction, and mid-solve consumers
        # (delay-equation history lookups especially) query near the current end
        forward = tspan[2] > tspan[1]
        tmid = tq[end ÷ 2]
        fresh = OrdinaryDiffEqCore.TsSearchHint(sol.t)
        @test fresh.idx_prev == 0
        @test OrdinaryDiffEqCore.ts_hint_start(fresh, sol.t) == lastindex(sol.t)
        @test OrdinaryDiffEqCore._searchsortedfirst(fresh, sol.t, tmid, 2, forward) ==
            OrdinaryDiffEqCore._searchsortedfirst(sol.t, tmid, 2, forward)
        @test fresh.idx_prev != 0
        @test OrdinaryDiffEqCore._searchsortedlast(fresh, sol.t, tmid, 1, forward) ==
            OrdinaryDiffEqCore._searchsortedlast(sol.t, tmid, 1, forward)
        # DDE-shaped: the grid grows in place after hint construction; a first
        # query near the new end is still exact
        tsg = sol.t[1:10]
        grown = OrdinaryDiffEqCore.TsSearchHint(tsg)
        append!(tsg, sol.t[11:end])
        @test OrdinaryDiffEqCore._searchsortedfirst(grown, tsg, tmid, 2, forward) ==
            OrdinaryDiffEqCore._searchsortedfirst(tsg, tmid, 2, forward)
        # correctness is strategy-independent: force each kind and re-check
        for kind in (
                OrdinaryDiffEqCore.KIND_BRACKET_GALLOP,
                OrdinaryDiffEqCore.KIND_INTERPOLATION_SEARCH,
            )
            hint.kind = kind
            @test all(sol(t) == ref[:left][t] for t in tq)
        end
    end
end

@testset "strategy selection from solve options" begin
    prob = ODEProblem(rot!, [1.0, 0.0], (0.0, 10.0))
    # fixed-dt stepping: uniform grid known at init, no re-probing
    sol = solve(prob, Tsit5(); adaptive = false, dt = 0.01)
    @test sol.interp.ts_hint.kind == OrdinaryDiffEqCore.KIND_INTERPOLATION_SEARCH
    @test sol.interp.ts_hint.probed_len == typemax(Int)
    # range saveat without save_everystep: uniform grid known at init
    sol = solve(prob, Tsit5(); saveat = 0.0:0.1:10.0)
    @test sol.interp.ts_hint.kind == OrdinaryDiffEqCore.KIND_INTERPOLATION_SEARCH
    # fully adaptive dense on a pure oscillator: dt settles into a narrow
    # controller band, so the ending-phase probe selects interpolation search
    osc!(du, u, p, t) = (du[1] = -u[2]; du[2] = u[1]; nothing)
    solu = solve(
        ODEProblem(osc!, [1.0, 0.0], (0.0, 1000.0)), Tsit5();
        abstol = 1.0e-8, reltol = 1.0e-8, dense = true
    )
    @test solu.interp.ts_hint.probed_len == length(solu.t)
    @test solu.interp.ts_hint.kind == OrdinaryDiffEqCore.KIND_INTERPOLATION_SEARCH
    # the damped rotation decays to steady state and dt grows ~25x: the grid
    # is bimodal, so the probe correctly keeps the gallop
    sold = solve(
        ODEProblem(rot!, [1.0, 0.0], (0.0, 1000.0)), Tsit5();
        abstol = 1.0e-8, reltol = 1.0e-8, dense = true
    )
    @test sold.interp.ts_hint.probed_len == length(sold.t)
    @test sold.interp.ts_hint.kind == OrdinaryDiffEqCore.KIND_BRACKET_GALLOP
end
