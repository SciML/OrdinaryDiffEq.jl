# The scalar `ode_interpolation`/`ode_interpolation!` paths warm-start their
# interval search from a `FindFirstFunctions.Guesser` stored in
# `InterpolationData`. These tests pin two properties:
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
        hint.idx_prev[] = 1
        sol(tq[end ÷ 2])
        @test hint.idx_prev[] != 1
    end
end
