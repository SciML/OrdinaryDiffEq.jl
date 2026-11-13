using GlobalDiffEq, OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, LinearAlgebra
# Loads the BrownFullBasicInit method.
import OrdinaryDiffEqNonlinearSolve
using DiffEqCallbacks: PresetTimeCallback, PeriodicCallback, SavingCallback, SavedValues,
    TerminateSteadyState
using Test
import SciMLBase, SciMLLogging

# DefectCorrection with callbacks, in both modes. The companion applies the user
# callbacks to the true state, so the estimate follows the true event times, which differ
# from the forward ones. Each case has an exact solution.
const T = 2.0
const MODES = (InterpolatingMode, SimultaneousMode)
decay!(du, u, p, t) = (du[1] = -u[1]; nothing)
decay(u, p, t) = [-u[1]]
switch!(du, u, p, t) = (du[1] = p[1] - u[1]; nothing)
decay_and_switch!(du, u, p, t) = (du[1] = -u[1]; du[2] = p[1] - u[2]; nothing)
# x' = -x, 0 = z - 2x: a mass-matrix DAE.
decay_dae!(du, u, p, t) = (du[1] = -u[1]; du[2] = u[2] - 2u[1]; nothing)
# x' = -x, 0 = z + z^3 - 2x: a mass-matrix DAE with a nonlinear constraint.
cubic_dae!(du, u, p, t) = (du[1] = -u[1]; du[2] = u[2] + u[2]^3 - 2u[1]; nothing)
function cubic_z(x)
    z = x
    for _ in 1:50
        z -= (z + z^3 - 2x) / (1 + 3z^2)
    end
    return z
end
reset_callback() = ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> (i.u[1] = 1.0; nothing))
# x' = -x from x = 1, with a reset to 1 at x = 0.5. The resets are at log 2 and 2 log 2.
const RESET_EXACT = exp(-(T - 2log(2)))

# A VectorContinuousCallback affect gets the index of the fired component, or a mask of
# the components that fired at the same time.
fired(index::Integer, k) = index == k
fired(mask, k) = !iszero(mask[k])
function reset_and_switch!(integrator, index)
    fired(index, 1) && (integrator.u[1] = 1.0)
    fired(index, 2) && (integrator.p[1] = 0.0)
    return nothing
end

# (name, problem, callback, exact u(T), extra solve keywords, forward solver)
const CASES = [
    ("reset", () -> ODEProblem(decay!, [1.0], (0.0, T)), reset_callback, [RESET_EXACT], (;), Tsit5()),
    (
        "reset, implicit solver", () -> ODEProblem(decay!, [1.0], (0.0, T)), reset_callback,
        [RESET_EXACT], (;), Rodas5P(),
    ),
    ("reset, out of place", () -> ODEProblem(decay, [1.0], (0.0, T)), reset_callback, [RESET_EXACT], (;), Tsit5()),
    (
        # The affect keeps the algebraic state consistent.
        "reset, mass-matrix DAE",
        () -> ODEProblem(
            ODEFunction(decay_dae!; mass_matrix = Diagonal([1.0, 0.0])), [1.0, 2.0], (0.0, T)
        ),
        () -> ContinuousCallback(
            (u, t, i) -> u[1] - 0.5, i -> (i.u[1] = 1.0; i.u[2] = 2.0; nothing)
        ),
        [RESET_EXACT, 2RESET_EXACT], (;), Rodas5P(),
    ),
    (
        # The affect resets only x, and the reinitialization solves the nonlinear
        # constraint for z. The companion must reinitialize the true state too.
        "reset, reinitialized nonlinear DAE",
        () -> ODEProblem(
            ODEFunction(cubic_dae!; mass_matrix = Diagonal([1.0, 0.0])),
            [1.0, cubic_z(1.0)], (0.0, T)
        ),
        () -> ContinuousCallback(
            (u, t, i) -> u[1] - 0.5, i -> (i.u[1] = 1.0; nothing);
            initializealg = BrownFullBasicInit()
        ),
        [RESET_EXACT, cubic_z(RESET_EXACT)], (;), Rodas5P(),
    ),
    (
        # The companion runs on a copy of the callback, so its counter starts at zero.
        "affect with internal state",
        () -> ODEProblem(decay!, [1.0], (0.0, T)),
        function ()
            fired = Ref(0)
            return ContinuousCallback(
                (u, t, i) -> u[1] - 0.5,
                i -> (fired[] == 0 && (i.u[1] = 1.0); fired[] += 1; nothing)
            )
        end,
        [exp(-(T - log(2)))], (;), Tsit5(),
    ),
    (
        # Up-crossing of x = 1 at log 2: jump to 1.5 and switch off the source.
        # Down-crossing at log 2 + log 1.5: switch the source to 0.5.
        "up- and down-crossing affects",
        () -> ODEProblem(switch!, [0.0], (0.0, T), [2.0]),
        () -> ContinuousCallback(
            (u, t, i) -> u[1] - 1.0,
            i -> (i.u[1] = 1.5; i.p[1] = 0.0; nothing), i -> (i.p[1] = 0.5; nothing)
        ),
        [0.5 + 0.5exp(-(T - log(2) - log(1.5)))], (;), Tsit5(),
    ),
    (
        "parameter switch", () -> ODEProblem(switch!, [0.0], (0.0, T), [2.0]),
        () -> ContinuousCallback((u, t, i) -> u[1] - 1.0, i -> (i.p[1] = 0.0; nothing)),
        [exp(-(T - log(2)))], (;), Tsit5(),
    ),
    (
        "vector callback", () -> ODEProblem(decay_and_switch!, [1.0, 0.0], (0.0, T), [2.0]),
        () -> VectorContinuousCallback(
            (out, u, t, i) -> (out[1] = u[1] - 0.5; out[2] = u[2] - 1.0; nothing),
            reset_and_switch!, 2
        ),
        [exp(-(T - 2log(2))), exp(-(T - log(2)))], (;), Tsit5(),
    ),
    (
        "fixed-time doubling", () -> ODEProblem(decay!, [1.0], (0.0, T)),
        () -> DiscreteCallback((u, t, i) -> t == 0.5, i -> (i.u[1] *= 2; nothing)),
        [2exp(-T)], (; tstops = [0.5]), Tsit5(),
    ),
    (
        "PresetTimeCallback", () -> ODEProblem(decay!, [1.0], (0.0, T)),
        () -> PresetTimeCallback(0.5, i -> (i.u[1] *= 2; nothing)), [2exp(-T)], (;), Tsit5(),
    ),
    (
        # Its state is in task-local storage. In SimultaneousMode, the copy in the companion
        # must not share that state with the forward solve. Events at 0.6, 1.2 and 1.8.
        "PeriodicCallback", () -> ODEProblem(decay!, [1.0], (0.0, T)),
        () -> PeriodicCallback(i -> (i.u[1] *= 1.5; nothing), 0.6), [1.5^3 * exp(-T)], (;), Tsit5(),
    ),
    (
        # The saved state at T is after the affect; P(T⁻) is before it.
        "event at the final time", () -> ODEProblem(decay!, [1.0], (0.0, T)),
        () -> DiscreteCallback((u, t, i) -> t == T, i -> (i.u[1] *= 2; nothing)),
        [2exp(-T)], (;), Tsit5(),
    ),
]

@testset "Estimate with events: $name, $mode" for (name, make_prob, make_callback, exact, extra, solver) in CASES,
        mode in MODES

    alg = GlobalErrorEstimation(solver; mode)
    for tol in (1.0e-4, 1.0e-6)
        prob = make_prob()
        p_before = deepcopy(prob.p)
        # In the VectorContinuousCallback case, both components fire at the same time. At
        # 1e-4 the forward solve can miss one of them, and the estimate then warns.
        est = global_error_estimate(
            prob, alg; callback = make_callback(), abstol = tol, reltol = tol, extra...
        )
        # The estimate does not change the user's parameters.
        @test prob.p == p_before
        sol = solve(make_prob(), solver; callback = make_callback(), abstol = tol, reltol = tol, extra...)
        @test est / norm(exact - sol.u[end]) ≈ 1 rtol = 1.0e-2
    end

    # The callback can also come from the problem keywords.
    prob = make_prob()
    in_problem = ODEProblem(prob.f, prob.u0, prob.tspan, prob.p; callback = make_callback())
    est = global_error_estimate(in_problem, alg; abstol = 1.0e-6, reltol = 1.0e-6, extra...)
    sol = solve(make_prob(), solver; callback = make_callback(), abstol = 1.0e-6, reltol = 1.0e-6, extra...)
    @test est / norm(exact - sol.u[end]) ≈ 1 rtol = 1.0e-2

    # Global error control with events.
    gtol = 1.0e-7
    sol = solve(
        make_prob(), GlobalErrorEstimation(solver; mode, gtol);
        callback = make_callback(), abstol = 1.0e-3, reltol = 1.0e-3, extra...
    )
    @test SciMLBase.successful_retcode(sol)
    @test norm(exact - sol.u[end]) <= gtol
end

@testset "A missed event is estimated and reported: $mode" for mode in MODES
    # Both components of the vector callback cross at t = log 2. At 1e-4 the forward solve
    # can fire only the first one, which depends on rounding and on the Julia version. Its
    # second state is then wrong by O(1); the companion fires both, and the estimator warns.
    prob = ODEProblem(decay_and_switch!, [1.0, 0.0], (0.0, T), [2.0])
    callback = VectorContinuousCallback(
        (out, u, t, i) -> (out[1] = u[1] - 0.5; out[2] = u[2] - 1.0; nothing),
        reset_and_switch!, 2
    )
    sol = solve(deepcopy(prob), Tsit5(); callback = deepcopy(callback), abstol = 1.0e-4, reltol = 1.0e-4)
    actual = norm([exp(-(T - 2log(2))), exp(-(T - log(2)))] - sol.u[end])
    alg = GlobalErrorEstimation(Tsit5(); mode)
    estimate() = global_error_estimate(prob, alg; callback, abstol = 1.0e-4, reltol = 1.0e-4)
    if actual > 0.1
        est = @test_logs (:warn, r"events of the forward solve") estimate()
        @test est / actual ≈ 1 rtol = 1.0e-3
    else
        est = @test_logs estimate()
        @test est / actual ≈ 1 rtol = 1.0e-2
    end
end

@testset "DiscreteCallback with a state condition: $mode" for mode in MODES
    # The condition is checked at the step ends of the forward solve, so the true
    # solution switches at the first forward step end where the true x is at least 1.
    x_before(t) = 2(1 - exp(-t))
    make_callback() = DiscreteCallback(
        (u, t, i) -> u[1] >= 1 && i.p[1] == 2.0, i -> (i.p[1] = 0.0; nothing)
    )
    alg = GlobalErrorEstimation(Tsit5(); mode)
    for tol in (1.0e-4, 1.0e-6)
        prob = ODEProblem(switch!, [0.0], (0.0, T), [2.0])
        est = global_error_estimate(prob, alg; callback = make_callback(), abstol = tol, reltol = tol)
        sol = solve(prob, Tsit5(); callback = make_callback(), abstol = tol, reltol = tol)
        step_ends = [unique(sol.t)[2:(end - 1)]; T]
        t_switch = step_ends[findfirst(t -> x_before(t) >= 1, step_ends)]
        actual = abs(x_before(t_switch) * exp(-(T - t_switch)) - sol.u[end][1])
        @test est / actual ≈ 1 rtol = 1.0e-2
    end
end

@testset "Error trajectory with several calls or one-sided saves: $mode" for mode in MODES
    # Each case doubles u at t = 0.5. A callback can save before its affect, after it, or
    # both; two callbacks at one time give four saves. Each entry is the state after the
    # calls before it, so its exact value follows from its number of doublings.
    double(i) = (i.u .*= 2; nothing)
    cases = (
        DiscreteCallback((u, t, i) -> t == 0.5, double; save_positions = (true, false)),
        DiscreteCallback((u, t, i) -> t == 0.5, double; save_positions = (false, true)),
        CallbackSet(
            DiscreteCallback((u, t, i) -> t == 0.5, double), DiscreteCallback((u, t, i) -> t == 0.5, double)
        ),
    )
    for callback in cases
        sol = solve(
            ODEProblem((u, p, t) -> -u, [1.0], (0.0, 1.0)), GlobalErrorEstimation(Tsit5(); mode);
            callback, tstops = [0.5], abstol = 1.0e-6, reltol = 1.0e-6
        )
        doublings = [round(Int, log2(u[1] / exp(-t))) for (t, u) in zip(sol.t, sol.u)]
        actual = [2.0^k * exp(-t) - u[1] for (k, t, u) in zip(doublings, sol.t, sol.u)]
        @test all(abs(e[1] - a) <= 1.0e-2 * abs(a) + 1.0e-14 for (e, a) in zip(sol.global_error, actual))
    end
end

@testset "Initialization affects: $mode" for mode in MODES
    # The forward initialization doubles u. The true solution starts with the same affect,
    # so the companion does not apply it to the true state again.
    callback = () -> PeriodicCallback(i -> (i.u .*= 2; nothing), 10.0; initial_affect = true)
    prob = ODEProblem(decay!, [1.0], (0.0, T))
    est = global_error_estimate(
        prob, GlobalErrorEstimation(Tsit5(); mode); callback = callback(), abstol = 1.0e-6, reltol = 1.0e-6
    )
    sol = solve(prob, Tsit5(); callback = callback(), abstol = 1.0e-6, reltol = 1.0e-6)
    @test est / abs(2exp(-T) - sol.u[end][1]) ≈ 1 rtol = 1.0e-2
end

@testset "Terminating callbacks: $mode" for mode in MODES
    # The solve ends at the forward event, and the estimate is the error at that time. The
    # true state must not terminate, even if the true event comes first.
    # TerminateSteadyState also reads the derivative of the true state.
    callbacks = (
        () -> ContinuousCallback((u, t, i) -> u[1] - 0.5, SciMLBase.terminate!),
        () -> TerminateSteadyState(1.0e-4, 1.0e-4),
    )
    # No warning: when the true event comes after the forward one, the companion does not
    # reach it, and that is not an event mismatch.
    for make_callback in callbacks, tol in (1.0e-4, 1.0e-6)
        prob = ODEProblem(decay!, [1.0], (0.0, 100.0))
        est = @test_logs min_level = Base.CoreLogging.Warn global_error_estimate(
            prob, GlobalErrorEstimation(Tsit5(); mode); callback = make_callback(), abstol = tol, reltol = tol
        )
        sol = solve(prob, Tsit5(); callback = make_callback(), abstol = tol, reltol = tol)
        @test sol.t[end] < 100
        @test est / abs(exp(-sol.t[end]) - sol.u[end][1]) ≈ 1 rtol = 1.0e-2
    end

    # At a forward step end, a DiscreteCallback terminates if x is below (or above) a value
    # between the forward and the true state. The two solutions decide differently, and a
    # warning reports it.
    tol = 1.0e-3
    prob = ODEProblem(decay!, [1.0], (0.0, 2.0))
    plain = solve(prob, Tsit5(); abstol = tol, reltol = tol)
    k = findfirst(>(0.7), plain.t)
    ts = plain.t[k]
    c = (plain.u[k][1] + exp(-ts)) / 2
    @test plain.u[k][1] > exp(-ts)
    for (below, message) in (
            (true, r"true solution.*terminates at t = .*but the forward solve does not"),
            (false, r"forward solve terminates at t = .*does not terminate at the same event"),
        )
        cb = DiscreteCallback(
            (u, t, i) -> t == ts,
            i -> ((below ? i.u[1] < c : i.u[1] > c) && SciMLBase.terminate!(i); nothing)
        )
        est = @test_logs (:warn, message) global_error_estimate(
            prob, GlobalErrorEstimation(Tsit5(); mode); callback = cb, abstol = tol, reltol = tol
        )
        # The `global_error_event_mismatch` toggle of `DEVerbosity` silences it.
        @test_logs min_level = Base.CoreLogging.Warn global_error_estimate(
            prob, GlobalErrorEstimation(Tsit5(); mode); callback = cb, abstol = tol, reltol = tol,
            verbose = DEVerbosity(global_error_event_mismatch = SciMLLogging.Silent())
        )
        sol = solve(prob, Tsit5(); callback = cb, abstol = tol, reltol = tol)
        @test sol.t[end] == (below ? 2.0 : ts)
        @test est / abs(exp(-sol.t[end]) - sol.u[end][1]) ≈ 1 rtol = 1.0e-2
    end
end

@testset "Unsupported callback functions fail clearly: $mode" for mode in MODES
    # An affect that changes the step size has no meaning for the true state.
    prob = ODEProblem(decay!, [1.0], (0.0, T))
    callback = ContinuousCallback(
        (u, t, i) -> u[1] - 0.5,
        i -> (i.u[1] = 1.0; SciMLBase.set_proposed_dt!(i, 0.01); nothing)
    )
    @test_throws ArgumentError global_error_estimate(
        prob, GlobalErrorEstimation(Tsit5(); mode); callback
    )
end

@testset "Companion reinitialization after an event" begin
    prob = ODEProblem(decay!, [1.0], (0.0, T))
    reinit(alg, kwargs = (;)) = GlobalDiffEq._companion_reinit(alg, prob, kwargs, 1.0e-12)
    # Without initialization data, the default reinitialization is CheckInit, which
    # changes nothing.
    @test reinit(nothing) === SciMLBase.NoInit()
    @test reinit(SciMLBase.CheckInit()) === SciMLBase.NoInit()
    # The initializealg of the solve applies when the callback has none. The companion uses
    # its own tolerance.
    @test reinit(BrownFullBasicInit(abstol = 1.0e-7)).abstol == 1.0e-12
    @test reinit(nothing, (; initializealg = BrownFullBasicInit())) isa BrownFullBasicInit
    # OverrideInit (the default for a ModelingToolkit problem) solves for u, not for the
    # companion state, so a DAE rejects it. An ODE does not reinitialize.
    @test_throws ArgumentError GlobalDiffEq._reinitialize_true_state!(
        (; isdae = true), SciMLBase.OverrideInit()
    )
    @test GlobalDiffEq._reinitialize_true_state!(
        (; isdae = false), SciMLBase.OverrideInit()
    ) === nothing
end

@testset "SimultaneousMode matches InterpolatingMode with events: $name" for (name, make_prob, make_callback, _, extra, solver) in CASES
    for tol in (1.0e-4, 1.0e-6)
        estimates = map(MODES) do mode
            global_error_estimate(
                make_prob(), GlobalErrorEstimation(solver; mode);
                callback = make_callback(), abstol = tol, reltol = tol, extra...
            )
        end
        @test estimates[2] ≈ estimates[1] rtol = 1.0e-3
    end
end

@testset "A DiscreteCallback that a continuous event triggers: $mode" for mode in MODES
    # At x = 0.5 the continuous event doubles x and sets a flag. At the same step end, a
    # DiscreteCallback adds x² to y and clears the flag. The true event ends a step of the
    # true solution, so the companion must check the DiscreteCallbacks there too.
    tied = () -> CallbackSet(
        ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> (i.u[1] *= 2; i.p[1] = 1.0; nothing)),
        DiscreteCallback((u, t, i) -> i.p[1] == 1.0, i -> (i.u[2] += i.u[1]^2; i.p[1] = 0.0; nothing))
    )
    prob = () -> ODEProblem((u, p, t) -> [-u[1], -u[2]], [1.0, 1.0], (0.0, 2.0), [0.0])
    ref = solve(prob(), Tsit5(); callback = tied(), abstol = 1.0e-13, reltol = 1.0e-13)
    for tol in (1.0e-5, 1.0e-7)
        sol = solve(prob(), Tsit5(); callback = tied(), abstol = tol, reltol = tol)
        est = global_error_estimate(
            prob(), GlobalErrorEstimation(Tsit5(); mode); callback = tied(), abstol = tol, reltol = tol
        )
        @test est / norm(ref.u[end] .- sol.u[end]) ≈ 1 rtol = 1.0e-3
    end
end

@testset "Callbacks at the forward event knots, integer time span, termination: $mode" for mode in MODES
    # The companion does not call a callback that acts at each step, such as a
    # SavingCallback, at the forward event time. The true solution ends a step at its own
    # event time instead. An extra call would give a false event-mismatch warning.
    reset = () -> ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> (i.u[1] = 1.0; nothing))
    saving = () -> SavingCallback((u, t, i) -> copy(u), SavedValues(Float64, Vector{Float64}))
    prob = ODEProblem((u, p, t) -> [-u[1]], [1.0], (0.0, 2.0))
    alg = GlobalErrorEstimation(Tsit5(); mode)
    est = @test_logs min_level = Base.CoreLogging.Warn global_error_estimate(
        prob, alg; callback = CallbackSet(reset(), saving()), abstol = 1.0e-6, reltol = 1.0e-6
    )
    # The copy of the SavingCallback in the companion acts at each knot, which can change
    # the companion steps slightly.
    est_reset = global_error_estimate(prob, alg; callback = reset(), abstol = 1.0e-6, reltol = 1.0e-6)
    @test est ≈ est_reset rtol = 1.0e-3

    # An integer time span gives the same estimate; the event times are floats.
    int_prob = ODEProblem((u, p, t) -> [-u[1]], [1.0], (0, 2))
    @test global_error_estimate(int_prob, alg; callback = reset(), abstol = 1.0e-6, reltol = 1.0e-6) ≈
        est_reset rtol = 1.0e-6

    # A terminating event triggers a DiscreteCallback at its step end. When the true event
    # comes after the forward end, the companion runs neither callback and does not warn.
    triggered = () -> CallbackSet(
        ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> (i.p[1] = 1.0; SciMLBase.terminate!(i))),
        DiscreteCallback((u, t, i) -> i.p[1] == 1.0, i -> (i.u[1] += 0.1; i.p[1] = 0.0; nothing))
    )
    flagged = () -> ODEProblem((u, p, t) -> [-u[1]], [1.0], (0.0, 2.0), [0.0])
    for tol in (1.0e-4, 1.0e-6)
        sol = solve(flagged(), Tsit5(); callback = triggered(), abstol = tol, reltol = tol)
        @test sol.t[end] < log(2)
        @test_logs min_level = Base.CoreLogging.Warn global_error_estimate(
            flagged(), alg; callback = triggered(), abstol = tol, reltol = tol
        )
    end
end

# A forward step a → b with b > 2a: a + (b - a) is not b, so the companion must stop at b
# itself, not at the end of a step of length b - a.
@testset "SimultaneousMode stops the companion at the forward step end" begin
    a, b = 0.008476769222961617, 0.34529383286546195
    @test a + (b - a) != b
    prob = ODEProblem((u, p, t) -> -0.1 .* u, [1.0, 2.0], (0.0, 1.0))
    kwargs = (; abstol = 1.0e-3, reltol = 1.0e-3, tstops = [a, b])
    @test b in solve(prob, Tsit5(); kwargs...).t
    estimates = map(MODES) do mode
        global_error_estimate(prob, GlobalErrorEstimation(Tsit5(); mode); kwargs...)
    end
    @test estimates[2] ≈ estimates[1] rtol = 1.0e-6
    trajectories = map(MODES) do mode
        solve(prob, GlobalErrorEstimation(Tsit5(); mode); kwargs...).global_error
    end
    @test trajectories[2] ≈ trajectories[1] rtol = 1.0e-6
end

@testset "SavingCallback sees only the forward solve: $mode, saveat = $saveat" for mode in MODES,
        saveat in ((), [0.3, 0.7, 1.0, 1.5, 1.9])

    saving(values) = CallbackSet(
        reset_callback(), SavingCallback((u, t, i) -> u[1], values; saveat)
    )
    make_prob() = ODEProblem(decay!, [1.0], (0.0, T))
    # The estimate saves the same values as a plain solve.
    values = SavedValues(Float64, Float64)
    est = global_error_estimate(
        make_prob(), GlobalErrorEstimation(Tsit5(); mode);
        callback = saving(values), abstol = 1.0e-6, reltol = 1.0e-6
    )
    plain = SavedValues(Float64, Float64)
    sol = solve(make_prob(), Tsit5(); callback = saving(plain), abstol = 1.0e-6, reltol = 1.0e-6)
    @test values.t == plain.t
    @test values.saveval == plain.saveval
    @test est / abs(RESET_EXACT - sol.u[end][1]) ≈ 1 rtol = 1.0e-2
    # Under global error control, only the final solve saves values. It repeats the
    # accepted solve, so it does not warn.
    values = SavedValues(Float64, Float64)
    sol = @test_logs min_level = Base.CoreLogging.Warn solve(
        make_prob(), GlobalErrorEstimation(Tsit5(); mode, gtol = 1.0e-7);
        callback = saving(values), abstol = 1.0e-3, reltol = 1.0e-3
    )
    @test SciMLBase.successful_retcode(sol)
    # The solution saves an event time twice (before and after the affect), the
    # SavingCallback once.
    @test values.t == (saveat == () ? unique(sol.t) : saveat)
end

@testset "A final forward solve that does not repeat the accepted one" begin
    # A callback with outside state, or a solve that is not deterministic, can make the
    # final forward solve differ from the accepted one. The test calls `_final_solve`
    # directly: a global counter in an affect would also move the true event, and the
    # control loop would not meet gtol.
    accepted = (; estimate = 0.0, sol = (; t = [0.0], u = [[1.0]]))
    values = SavedValues(Float64, Float64)
    callback = CallbackSet(reset_callback(), SavingCallback((u, t, i) -> u[1], values))
    sol = @test_logs (:warn, r"does not repeat the accepted one") GlobalDiffEq._final_solve(
        ODEProblem(decay!, [1.0], (0.0, T)), GlobalErrorEstimation(Tsit5()), accepted;
        callback, abstol = 1.0e-6, reltol = 1.0e-6
    )
    @test length(sol.global_error) == length(sol.t)
    # The new estimate runs on a copy of the callbacks made before the forward solve, so
    # the SavingCallback of the user has the values of the forward solve only.
    @test values.t == unique(sol.t)
    # The `global_error_final_solve` toggle of `DEVerbosity` silences the warning.
    @test_logs min_level = Base.CoreLogging.Warn GlobalDiffEq._final_solve(
        ODEProblem(decay!, [1.0], (0.0, T)), GlobalErrorEstimation(Tsit5()), accepted;
        callback = reset_callback(), abstol = 1.0e-6, reltol = 1.0e-6,
        verbose = DEVerbosity(global_error_final_solve = SciMLLogging.Silent())
    )
    @test sol.global_error[end][1] / (RESET_EXACT - sol.u[end][1]) ≈ 1 rtol = 1.0e-2
end

@testset "Error trajectory with events: $mode" for mode in MODES
    # The true solution resets at log 2 and 2 log 2, and the forward solution at nearby
    # times. The solution saves each forward event time twice. Between a true and a forward
    # event, the error is 0.5.
    reset_x(t) = exp(-(t - fld(t, log(2)) * log(2)))
    for tol in (1.0e-4, 1.0e-6)
        sol = solve(
            ODEProblem(decay!, [1.0], (0.0, T)), GlobalErrorEstimation(Tsit5(); mode);
            callback = reset_callback(), abstol = tol, reltol = tol
        )
        actual = [reset_x(t) - u[1] for (t, u) in zip(sol.t, sol.u)]
        @test maximum(abs, actual) ≈ 0.5 rtol = 1.0e-3
        @test all(abs(e[1] - a) <= 1.0e-3 * tol for (e, a) in zip(sol.global_error, actual))
    end

    # A fixed-time doubling: the true solution jumps at the same time, so the entry before
    # the event and the entry after it have different errors.
    sol = solve(
        ODEProblem(decay!, [1.0], (0.0, T)), GlobalErrorEstimation(Tsit5(); mode);
        callback = PresetTimeCallback([1.0], i -> (i.u[1] *= 2; nothing)),
        abstol = 1.0e-6, reltol = 1.0e-6
    )
    before_doubling(i) = sol.t[i] < 1 || (sol.t[i] == 1 && sol.t[i + 1] == 1)
    actual = [(before_doubling(i) ? 1 : 2) * exp(-sol.t[i]) - sol.u[i][1] for i in eachindex(sol.t)]
    @test count(==(1.0), sol.t) == 2
    @test all(
        abs(e[1] - a) <= 1.0e-2 * maximum(abs, actual) for (e, a) in zip(sol.global_error, actual)
    )

    # A DAE with a reinitialization after the event: both states have an error trajectory.
    sol = solve(
        ODEProblem(
            ODEFunction(cubic_dae!; mass_matrix = Diagonal([1.0, 0.0])),
            [1.0, cubic_z(1.0)], (0.0, T)
        ),
        GlobalErrorEstimation(Rodas5P(); mode);
        callback = ContinuousCallback(
            (u, t, i) -> u[1] - 0.5, i -> (i.u[1] = 1.0; nothing);
            initializealg = BrownFullBasicInit()
        ),
        abstol = 1.0e-6, reltol = 1.0e-6
    )
    actual = [[reset_x(t), cubic_z(reset_x(t))] .- u for (t, u) in zip(sol.t, sol.u)]
    @test maximum(norm, actual) > 0.1
    @test all(norm(e - a) <= 1.0e-6 for (e, a) in zip(sol.global_error, actual))
end
