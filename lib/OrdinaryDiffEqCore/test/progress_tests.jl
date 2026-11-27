using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqCore: progress_fraction, progress_min_interval_ns, PROGRESS_TSPAN_RTOL
using SciMLBase, Logging, Test

progress_records(logger) = [r for r in logger.logs if haskey(r.kwargs, :progress)]
intermediate_fracs(logger) = [
    r.kwargs[:progress] for r in progress_records(logger)
        if r.kwargs[:progress] isa Real && !iszero(r.kwargs[:progress])
]

const THROTTLE_OPEN = UInt64(0)
const THROTTLE_CLOSED = typemax(UInt64)

# Solve with a fixed wall-clock throttle, so that the outcome never depends on how fast
# the machine is: with `THROTTLE_OPEN` only the tspan-fraction gate decides, and with
# `THROTTLE_CLOSED` no intermediate message can get through.
function progress_solve(prob, alg; throttle, kwargs...)
    logger = Test.TestLogger(min_level = LogLevel(-1))
    integrator = with_logger(logger) do
        integrator = init(prob, alg; progress = true, kwargs...)
        integrator.opts.progress_min_interval_ns = throttle
        solve!(integrator)
        integrator
    end
    return integrator, logger
end

function robertson(du, u, p, t)
    y1, y2, y3 = u
    k1, k2, k3 = p
    du[1] = -k1 * y1 + k3 * y2 * y3
    du[2] = k1 * y1 - k2 * y2^2 - k3 * y2 * y3
    du[3] = k2 * y2^2
    return nothing
end
rober_prob = ODEProblem(robertson, [1.0, 0.0, 0.0], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))

@testset "progress_fraction" begin
    @test progress_fraction(0.0, (0.0, 10.0)) == 0.0
    @test progress_fraction(5.0, (0.0, 10.0)) == 0.5
    @test progress_fraction(10.0, (0.0, 10.0)) == 1.0
    @test progress_fraction(15.0, (5.0, 25.0)) == 0.5
    # reverse time
    @test progress_fraction(7.5, (10.0, 0.0)) == 0.25
    # outside tspan (reinit!, add_tstop! past the end): no fraction
    @test isnan(progress_fraction(11.0, (0.0, 10.0)))
    @test isnan(progress_fraction(-1.0, (0.0, 10.0)))
    @test isnan(progress_fraction(-1.0, (10.0, 0.0)))
    @test isnan(progress_fraction(NaN, (0.0, 10.0)))
    # not computed there, where the time type may not represent it (10^20 here)
    @test isnan(progress_fraction(10^10 // 1, (0 // 1, 1 // 10^10)))
    @test progress_fraction(1 // (2 * 10^10), (0 // 1, 1 // 10^10)) == 0.5
    # undefined fraction
    @test isnan(progress_fraction(1.0, (0.0, Inf)))
    @test isnan(progress_fraction(0.0, (0.0, 0.0)))
    @test isnan(progress_fraction(0 // 1, (0 // 1, 0 // 1)))  # no DivideError
    @test progress_fraction(1.0f0, (0.0f0, 2.0f0)) isa Float64
    @test progress_fraction(1, (0, 4)) === 0.25
    # computed in the time type, so large offsets keep their precision
    @test progress_fraction(big"1e20" + 0.5, (big"1e20", big"1e20" + 1)) == 0.5
    # and spans outside the Float64 range work
    @test progress_fraction(big"5e-401", (big"0", big"1e-400")) ≈ 0.5
    @test progress_fraction(big"5e399", (big"0", big"1e400")) ≈ 0.5
    # a time type that can't be converted to Float64 gives NaN (no updates), not an error
    struct NoFloat64 <: AbstractFloat end
    @test isnan(OrdinaryDiffEqCore._progress_float(NoFloat64()))
end

@testset "Default throttle" begin
    @test progress_min_interval_ns() == 20_000_000
    integrator = init(ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0)), Tsit5(); progress = true)
    @test integrator.opts.progress_min_interval_ns == progress_min_interval_ns()
    @test integrator.opts.progress_tspan == OrdinaryDiffEqCore.DEFAULT_PROGRESS_TSPAN
    @test integrator.opts.progress_steps == 1
    @test init(ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0)), Tsit5()).opts.progress_min_interval_ns ==
        progress_min_interval_ns()
    # the throttle clock starts at init, and only with progress on
    @test integrator.opts.progress_last_ns > 0
    @test init(ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0)), Tsit5()).opts.progress_last_ns == 0
end

@testset "Throttle holds between messages" begin
    # Open the throttle exactly once: with an interval of half the UInt64 range, a last
    # message "one interval ago" (wrapping arithmetic) lets the first update through, and
    # every later one falls inside the interval however long the solve takes.
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 10.0))
    logger = Test.TestLogger(min_level = LogLevel(-1))
    integrator = with_logger(logger) do
        integrator = init(prob, Tsit5(); progress = true)
        interval = typemax(UInt64) ÷ 2
        integrator.opts.progress_min_interval_ns = interval
        integrator.opts.progress_last_ns = time_ns() - interval
        solve!(integrator)
        integrator
    end
    @test length(intermediate_fracs(logger)) == 1
    @test integrator.opts.progress_last_frac >= 1 - integrator.opts.progress_tspan
end

@testset "The throttle clock starts at init" begin
    # Not at the first accepted step: when that step already passes the tspan gate
    # (here 10% of tspan), its update must not be dropped however long the step took.
    # The sleep stands in for a slow first step; it guarantees at least 2 ms on the
    # clock, so with a 1 ms throttle the outcome doesn't depend on machine speed.
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    logger = Test.TestLogger(min_level = LogLevel(-1))
    integrator = with_logger(logger) do
        integrator = init(prob, Tsit5(); progress = true, adaptive = false, dt = 0.1)
        integrator.opts.progress_min_interval_ns = UInt64(1_000_000)
        sleep(0.002)
        step!(integrator)
        integrator
    end
    @test intermediate_fracs(logger) ≈ [0.1]

    # reinit! restarts the clock
    integrator.opts.progress_last_ns = UInt64(1)
    reinit!(integrator)
    @test integrator.opts.progress_last_ns > 1
end

@testset "Rejected attempts never emit" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 10.0))
    logger = Test.TestLogger(min_level = LogLevel(-1))
    with_logger(logger) do
        integrator = init(prob, Tsit5(); progress = true)
        integrator.opts.progress_min_interval_ns = THROTTLE_OPEN
        step!(integrator)
        empty!(logger.logs)  # the step itself may have sent an update
        integrator.opts.progress_last_frac = 0.0  # as if the gate was skipped (progress_steps > 1)
        integrator.accept_step = false
        OrdinaryDiffEqCore.maybe_log_progress!(integrator)
        @test isempty(intermediate_fracs(logger))
        integrator.accept_step = true
        OrdinaryDiffEqCore.maybe_log_progress!(integrator)
        @test length(intermediate_fracs(logger)) == 1
    end
end

@testset "DEOptions constructors for the previous field list" begin
    opts = init(ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0)), Tsit5()).opts
    new_fields = (
        :progress_tspan, :progress_last_frac, :progress_last_ns, :progress_min_interval_ns,
    )
    old_args = [getfield(opts, f) for f in fieldnames(typeof(opts)) if f ∉ new_fields]
    @test length(old_args) == 43
    for o in (typeof(opts)(old_args...), OrdinaryDiffEqCore.DEOptions(old_args...))
        @test o isa typeof(opts)
        @test o.progress_tspan == OrdinaryDiffEqCore.DEFAULT_PROGRESS_TSPAN
        @test o.progress_min_interval_ns == progress_min_interval_ns()
        @test o.progress_last_ns == 0
        @test o.progress_id == opts.progress_id
        @test o.timeseries_errors == opts.timeseries_errors
    end
end

@testset "Closed throttle: only the 0 and done messages, gate still re-arms" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    integrator, logger = progress_solve(prob, Tsit5(); throttle = THROTTLE_CLOSED)
    recs = progress_records(logger)
    @test first(recs).kwargs[:progress] == 0
    @test last(recs).kwargs[:progress] == "done"
    @test length(recs) == 2
    # the tspan gate kept passing even though every message was dropped
    @test integrator.opts.progress_last_frac >= 1 - integrator.opts.progress_tspan
end

@testset "progress_tspan = 0 leaves only the throttle" begin
    # one message per accepted step; rejected steps don't advance t and never emit
    integrator, logger = progress_solve(
        rober_prob, Rodas5P(); throttle = THROTTLE_OPEN, progress_tspan = 0
    )
    @test integrator.stats.nreject > 0
    @test length(intermediate_fracs(logger)) == integrator.stats.naccept
end

@testset "tspan-fraction gate" begin
    for (prob, alg) in (
                (ODEProblem((u, p, t) -> -u, 1.0, (0.0, 10.0)), Tsit5()),
                (ODEProblem((u, p, t) -> -u, 1.0, (10.0, 0.0)), Tsit5()),  # reverse time
                (ODEProblem((u, p, t) -> -u, 1.0f0, (0.0f0, 10.0f0)), Tsit5()),  # Float32 time
                (rober_prob, Rodas5P()),
            ),
            progress_tspan in (0.01, 0.1)

        _, logger = progress_solve(prob, alg; throttle = THROTTLE_OPEN, progress_tspan)
        fracs = intermediate_fracs(logger)
        @test !isempty(fracs)
        @test all(f -> f isa Float64, fracs)  # what ProgressLogging accepts
        @test all(0 .<= fracs .<= 1)
        @test issorted(fracs)
        @test all(diff(fracs) .>= progress_tspan * (1 - PROGRESS_TSPAN_RTOL))
        @test first(fracs) >= progress_tspan * (1 - PROGRESS_TSPAN_RTOL)
        @test length(fracs) <= 1 / progress_tspan
        @test last(progress_records(logger)).kwargs[:progress] == "done"
    end
end

@testset "Steps that divide progress_tspan are not skipped" begin
    # fixed dt of exactly 1% of tspan: one update per step, not every other step, also
    # with the rounding of Float32 time
    for T in (Float64, Float32)
        integrator, logger = progress_solve(
            ODEProblem((u, p, t) -> -u, one(T), (zero(T), one(T))), Tsit5();
            throttle = THROTTLE_OPEN, adaptive = false, dt = T(0.01)
        )
        @test length(intermediate_fracs(logger)) >= 99
    end
end

@testset "Robertson: updates spread over tspan, not bunched early" begin
    _, logger = progress_solve(rober_prob, Rodas5P(); throttle = THROTTLE_OPEN)
    fracs = intermediate_fracs(logger)
    # Rodas5P takes 60 steps here, 60% of them in the first 1% of tspan, so the old
    # every-1000-steps trigger sent no intermediate update; now they cover the span
    @test maximum(fracs) > 0.9
    @test length(fracs) >= 10
end

@testset "progress_steps still gates" begin
    _, logger = progress_solve(
        rober_prob, Rodas5P(); throttle = THROTTLE_OPEN, progress_steps = typemax(Int)
    )
    @test isempty(intermediate_fracs(logger))
end

@testset "reinit! resets the progress state" begin
    integrator, _ = progress_solve(
        ODEProblem((u, p, t) -> -u, 1.0, (0.0, 10.0)), Tsit5(); throttle = THROTTLE_OPEN
    )
    @test integrator.opts.progress_last_frac > 0
    last_ns = integrator.opts.progress_last_ns
    reinit!(integrator)
    @test integrator.opts.progress_last_frac == 0
    @test integrator.opts.progress_last_ns >= last_ns  # clock restarted

    # Progress stays measured against prob.tspan after reinit! onto another interval,
    # and only the part of the interval inside prob.tspan sends updates (no fractions
    # below 0 or above 1). Outside it, the gate restarts from 0.
    integrator = init(
        ODEProblem((u, p, t) -> -u, 1.0, (0.0, 10.0)), Tsit5();
        progress = true, adaptive = false, dt = 1.0e-3
    )
    integrator.opts.progress_min_interval_ns = THROTTLE_OPEN
    for (t0, tf, nmin, nmax) in (
            (10.0, 20.0, 0, 0), (-20.0, -10.0, 0, 0), (-5.0, 5.0, 45, 50), (5.0, 15.0, 45, 50),
        )
        reinit!(integrator; t0, tf)
        @test integrator.opts.progress_last_frac ==
            (0 <= t0 <= 10 ? progress_fraction(t0, (0.0, 10.0)) : 0.0)
        logger = Test.TestLogger(min_level = LogLevel(-1))
        with_logger(() -> solve!(integrator), logger)
        fracs = intermediate_fracs(logger)
        @test nmin <= length(fracs) <= nmax
        @test all(0 .< fracs .<= 1)
        @test all(diff(fracs) .>= integrator.opts.progress_tspan * (1 - PROGRESS_TSPAN_RTOL))
        @test last(progress_records(logger)).kwargs[:progress] == "done"
    end

    # an empty rational tspan must not break reinit!, with or without progress
    for progress in (false, true)
        integrator = init(
            ODEProblem((u, p, t) -> -u, 1.0, (0 // 1, 0 // 1)), Tsit5();
            dt = 1 // 10, adaptive = false, progress
        )
        reinit!(integrator)  # reached without a DivideError
        @test integrator.opts.progress_last_frac == 0
    end

    # reinit! far outside a Rational tspan, where the fraction would overflow Int64
    # (Tsit5 can't step in Rational time, so the per-step check is called directly)
    for progress in (false, true)
        integrator = init(
            ODEProblem((u, p, t) -> zero(u), 1.0, (0 // 1, 1 // 10^10)), Tsit5();
            adaptive = false, dt = 1 // 1, progress
        )
        integrator.opts.progress_min_interval_ns = THROTTLE_OPEN
        reinit!(integrator; t0 = 10^10 // 1, tf = (10^10 + 1) // 1)
        @test integrator.t == 10^10
        @test integrator.opts.progress_last_frac == 0
        logger = Test.TestLogger(min_level = LogLevel(-1))
        with_logger(() -> OrdinaryDiffEqCore.maybe_log_progress!(integrator), logger)
        @test isempty(progress_records(logger))
    end
end

@testset "No updates past the end of tspan" begin
    # a stop added beyond tspan[2]: the solve goes there, but sends nothing above 1
    added = Ref(false)
    cb = DiscreteCallback(
        (u, t, integrator) -> !added[],
        integrator -> (added[] = true; add_tstop!(integrator, 1.25))
    )
    integrator, logger = progress_solve(
        ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0)), Tsit5();
        throttle = THROTTLE_OPEN, adaptive = false, dt = 0.01, callback = cb
    )
    @test integrator.t == 1.25
    fracs = intermediate_fracs(logger)
    @test 95 <= length(fracs) <= 100
    @test all(0 .< fracs .<= 1)
    @test last(progress_records(logger)).kwargs[:progress] == "done"
end

@testset "progress_tspan and progress_min_interval are validated" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    for bad in (-0.1, NaN, Inf)
        @test_throws ArgumentError init(prob, Tsit5(); progress = true, progress_tspan = bad)
    end
    @test init(prob, Tsit5(); progress_tspan = 0).opts.progress_tspan === 0.0
    for bad in (-0.1, NaN)
        @test_throws ArgumentError init(prob, Tsit5(); progress = true, progress_min_interval = bad)
    end
    interval(x) = init(prob, Tsit5(); progress_min_interval = x).opts.progress_min_interval_ns
    @test interval(0) == 0
    @test interval(0.05) == 50_000_000
    @test interval(Inf) == typemax(UInt64)
end

@testset "progress_min_interval = 0 restores every progress_steps step" begin
    logger = Test.TestLogger(min_level = LogLevel(-1))
    integrator = with_logger(logger) do
        solve!(
            init(
                ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0)), Tsit5();
                progress = true, progress_tspan = 0, progress_min_interval = 0,
                adaptive = false, dt = 0.001, progress_steps = 10
            )
        )
    end
    @test length(intermediate_fracs(logger)) == 100
end

@testset "Infinite tspan: only the 0 and done messages" begin
    # the fraction is undefined, so no intermediate update is sent even with the throttle open
    stop = DiscreteCallback((u, t, integrator) -> t >= 5, terminate!)
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, Inf))
    _, logger = progress_solve(
        prob, Tsit5(); throttle = THROTTLE_OPEN, progress_tspan = 0,
        adaptive = false, dt = 0.01, callback = stop
    )
    @test [r.kwargs[:progress] for r in progress_records(logger)] == [0, "done"]
end
