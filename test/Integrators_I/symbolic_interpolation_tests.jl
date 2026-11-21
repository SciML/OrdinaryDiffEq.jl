using OrdinaryDiffEq, Test
using SciMLBase: SymbolCache

# `SymbolCache` resolves an `Expr` as an observed quantity, so `:(x + y)` is not a
# component of the state vector and cannot be produced by the dense-output interpolant.
f = ODEFunction(
    (du, u, p, t) -> (du[1] = p[1] * u[2]; du[2] = -p[1] * u[1]);
    sys = SymbolCache([:x, :y], [:ω], :t)
)
prob = ODEProblem(f, [1.0, 0.0], (0.0, 10.0), [2.0])

@testset "Symbolic idxs in the current-step interpolant" begin
    integrator = init(prob, Tsit5())
    step!(integrator, 1.0, true)
    mid = (integrator.tprev + integrator.t) / 2

    # A symbolic index evaluates the whole-state interpolant and then projects, while an
    # integer index evaluates the interpolant for that component alone. The two orders of
    # operations agree to rounding rather than bitwise, so these compare with `≈`.
    @test integrator(mid; idxs = :x) ≈ integrator(mid; idxs = 1)
    @test integrator(mid; idxs = :y) ≈ integrator(mid; idxs = 2)
    @test integrator(mid; idxs = [:x, :y]) ≈ integrator(mid)
    @test integrator(mid; idxs = :(x + y)) ≈ sum(integrator(mid))

    ts = [integrator.tprev, mid, integrator.t]
    @test integrator(ts; idxs = :x).u ≈ [integrator(tt; idxs = 1) for tt in ts]
    @test integrator(ts; idxs = :(x + y)).u ≈ [sum(integrator(tt)) for tt in ts]
    @test integrator(ts; idxs = [:x, :(x + y)]).t == ts

    # Symbolic routing does not disturb the existing integer and `nothing` paths
    @test integrator(mid) == integrator(mid; idxs = nothing)
    @test integrator(mid; idxs = [1, 2]) ≈ integrator(mid)
    @test integrator(ts; idxs = 1) == [integrator(tt; idxs = 1) for tt in ts]

    @test integrator(mid, Val{1}; idxs = :x) ≈ integrator(mid, Val{1}; idxs = 1)
    @test_throws ErrorException integrator(mid, Val{1}; idxs = :(x + y))
end
