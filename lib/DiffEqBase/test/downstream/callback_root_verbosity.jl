using OrdinaryDiffEq, SciMLBase, Test
using DiffEqBase: DEVerbosity, root_verbosity
using SciMLLogging: None, Standard

# u' = 1 from u(t0) = t0 gives u(t) = t, which Tsit5 and its interpolant reproduce exactly, so the
# event of u - c is at t = c up to rounding.
f(u, p, t) = 1.0
const VERBOSITIES = (DEVerbosity(None()), DEVerbosity(Standard()), DEVerbosity(), nothing)
kw(v) = v === nothing ? (;) : (; verbose = v)

@testset "root_verbosity of a real integrator" begin
    for v in VERBOSITIES
        integ = init(ODEProblem(f, 0.0, (0.0, 1.0)), Tsit5(); kw(v)...)
        expected = something(v, DEVerbosity()).nonlinear_verbosity
        @test root_verbosity(integ) == (; verbose = expected)
        @test (@inferred root_verbosity(integ)) isa NamedTuple{(:verbose,), Tuple{typeof(expected)}}
    end
end

@testset "event times with the integrator's verbosity" begin
    for v in VERBOSITIES
        # scalar, forward
        tev = Float64[]
        cb = ContinuousCallback((u, t, i) -> u - 0.3, i -> push!(tev, i.t))
        solve(ODEProblem(f, 0.0, (0.0, 1.0)), Tsit5(); callback = cb, kw(v)...)
        @test length(tev) == 1 && abs(tev[1] - 0.3) <= 1.0e-12

        # vector, forward and reverse time
        for (tspan, order) in (((0.0, 1.0), [1, 2]), ((1.0, 0.0), [2, 1]))
            hits = Tuple{Int, Float64}[]
            vcb = VectorContinuousCallback(
                (out, u, t, i) -> (out[1] = u - 0.3; out[2] = u - 0.8),
                (i, events) -> foreach(((idx, d),) -> iszero(d) || push!(hits, (idx, i.t)), enumerate(events)), 2
            )
            solve(ODEProblem(f, tspan[1], tspan), Tsit5(); callback = vcb, kw(v)...)
            @test first.(hits) == order
            @test all(abs(t - (idx == 1 ? 0.3 : 0.8)) <= 1.0e-12 for (idx, t) in hits)
        end

        # successive events in one step, which take the repeated-event bracket branch
        record = Int[]
        cb0 = ContinuousCallback((u, t, i) -> t - 0.0, i -> push!(record, 0); abstol = 0.0)
        vcb = VectorContinuousCallback(
            (out, u, t, i) -> out .= (t - 1.0e-8, t - 2.0e-8, t - 2.0e-7),
            (i, events) -> foreach(((idx, d),) -> iszero(d) || push!(record, idx), enumerate(events)),
            3; abstol = 0.0
        )
        solve(ODEProblem(f, 0.0, (-1.0, 1.0)), Tsit5(); dt = 2.0, callback = CallbackSet(cb0, vcb), kw(v)...)
        @test record == [0, 1, 2, 3]
    end
end
