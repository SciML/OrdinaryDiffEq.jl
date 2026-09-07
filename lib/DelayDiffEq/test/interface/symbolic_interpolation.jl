using DelayDiffEq
using OrdinaryDiffEqTsit5
using SciMLBase: SymbolCache
using Test

# `SymbolCache` resolves an `Expr` as an observed quantity, so `:(x + y)` is not a component
# of the state vector and cannot be produced by the dense-output interpolant.
f = DDEFunction(
    (du, u, h, p, t) -> (du[1] = -h(p, t - 0.2)[1] + u[2]; du[2] = -u[1]);
    sys = SymbolCache([:x, :y], [:τ], :t)
)
prob = DDEProblem(f, ones(2), (p, t) -> zeros(2), (0.0, 10.0), [0.2]; constant_lags = [0.2])

@testset "Symbolic idxs in the current-step interpolant" begin
    integrator = init(prob, MethodOfSteps(Tsit5()))
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

    # Symbolic routing does not disturb the existing integer and `nothing` paths
    @test integrator(mid) == integrator(mid; idxs = nothing)
    @test integrator(mid; idxs = [1, 2]) ≈ integrator(mid)

    @test_throws ErrorException integrator(mid, Val{1}; idxs = :(x + y))
end


# `HistoryODEIntegrator` is deliberately not routed: it implements none of the
# `SymbolicIndexingInterface` methods a `DEIntegrator` needs (it has no `p` field), so
# symbolic indexing against it is unsupported here and on master alike.
