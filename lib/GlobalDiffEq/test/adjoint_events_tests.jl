using GlobalDiffEq, OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using SciMLSensitivity, LinearAlgebra, Test
using DiffEqCallbacks: PresetTimeCallback, PeriodicCallback, SavingCallback, SavedValues,
    TerminateSteadyState
import ADTypes

relerr(estimate, actual) = abs(estimate - actual) / abs(actual)

# The estimate with one projection per state (the 2-norm of the linearized estimate, with
# no sampling error) and the actual error of a plain solve at its final time (on a copy: an
# affect can change the parameters of the problem in place).
function compare(prob, alg, tol, exact; adjoint_alg = Tsit5())
    ga = GlobalAdjoint(alg; samples = length(prob.u0), adjoint_alg)
    estimate = adjoint_error_estimate(prob, ga; abstol = tol, reltol = tol)
    sol = solve(deepcopy(prob), alg; abstol = tol, reltol = tol)
    return estimate, norm(exact(sol.t[end]) .- sol.u[end])
end

@testset "GlobalAdjoint: jumps of the dense output" begin
    # Rosenbrock32's dense output does not end at the step's state, so the identity needs
    # the knot terms and the endpoint shift.
    decay = ODEProblem((u, p, t) -> -u, [1.0], (0.0, 3.0))
    @test relerr(compare(decay, Rosenbrock32(), 1.0e-5, t -> [exp(-t)])...) < 1.0e-5
end

@testset "GlobalAdjoint: events" begin
    lv = ODEProblem(
        (u, p, t) -> [1.5u[1] - u[1] * u[2], -3.0u[2] + u[1] * u[2]], [1.0, 1.0], (0.0, 10.0);
        callback = ContinuousCallback((u, t, i) -> u[1] - 3.0, i -> (i.u[1] = 1.5; nothing), nothing)
    )
    lvref = solve(lv, Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    # A linearization: the relative error falls with the tolerance.
    @test relerr(compare(lv, Tsit5(), 1.0e-3, lvref)...) < 1.0e-2
    @test relerr(compare(lv, Tsit5(), 1.0e-5, lvref)...) < 1.0e-3

    period = log(55)
    slow = ODEProblem(
        (u, p, t) -> [0.45 - u[1]], [1.0], (0.0, 10.0);
        callback = ContinuousCallback((u, t, i) -> u[1] - 0.46, nothing, i -> (i.u[1] = 1.0; nothing))
    )
    slow_exact = t -> [0.45 + 0.55exp(-(t - fld(t, period) * period))]
    @test relerr(compare(slow, Tsit5(), 1.0e-3, slow_exact)...) < 1.0e-3

    preset = ODEProblem(
        (u, p, t) -> -u, [1.0], (0.0, 3.5);
        callback = DiscreteCallback((u, t, i) -> t in (1.0, 2.0, 3.0), i -> (i.u .*= 2; nothing)),
        tstops = [1.0, 2.0, 3.0]
    )
    @test relerr(compare(preset, Tsit5(), 1.0e-4, t -> [exp(-t) * 2.0^min(3, floor(Int, t))])...) < 1.0e-4

    walls = VectorContinuousCallback(
        (out, u, t, i) -> (out[1] = u[1] - 0.8; out[2] = u[1] + 0.8; nothing),
        (i, idx) -> (i.u[2] = -0.9 * i.u[2]; nothing), 2
    )
    wallprob = ODEProblem((u, p, t) -> [u[2], -u[1]], [0.0, 1.0], (0.0, 10.0); callback = walls)
    wallref = solve(wallprob, Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    @test relerr(compare(wallprob, Tsit5(), 1.0e-5, wallref)...) < 1.0e-3

    # Global error control with events.
    controlled = solve(lv, GlobalAdjoint(Tsit5(); gtol = 1.0e-5, samples = 2, adjoint_alg = Tsit5()))
    @test norm(lvref(10.0) .- controlled.u[end]) < 2.0e-5
end

@testset "GlobalAdjoint: DAEs" begin
    M = Diagonal([1.0, 0.0])
    fd = Rodas5P(autodiff = ADTypes.AutoFiniteDiff())
    cubic_z(x) = (q = sqrt(x^2 / 4 + 1 / 27); cbrt(x / 2 + q) + cbrt(x / 2 - q))
    cubic!(du, u, p, t) = (du[1] = -u[1]; du[2] = u[2]^3 + u[2] - u[1]; nothing)
    smooth = ODEProblem(ODEFunction(cubic!; mass_matrix = M), [2.0, 1.0], (0.0, 5.0))
    estimate, actual = compare(smooth, Rodas5P(), 1.0e-4, t -> (x = 2exp(-t); [x, cubic_z(x)]); adjoint_alg = fd)
    @test relerr(estimate, actual) < 1.0e-4
    per = log(4.0)
    reset = ODEProblem(
        ODEFunction(cubic!; mass_matrix = M), [2.0, 1.0], (0.0, 5.0);
        callback = ContinuousCallback(
            (u, t, i) -> u[1] - 0.5, nothing, i -> (i.u[1] = 2.0; nothing);
            initializealg = BrownFullBasicInit()
        )
    )
    exact = t -> (x = 2exp(-(t - fld(t, per) * per)); [x, cubic_z(x)])
    @test relerr(compare(reset, Rodas5P(), 1.0e-4, exact; adjoint_alg = fd)...) < 1.0e-3
end

@testset "GlobalAdjoint: events that change the parameters" begin
    # The rate changes from 1 to 3 when u reaches 1/2. The dual uses the parameters of each
    # interval, and the saltation matrix the vector field after the event (S = φ⁺/φ⁻ = 3).
    ts = log(2.0)
    switch_exact = t -> [t < ts ? exp(-t) : 0.5exp(-3(t - ts))]
    switch = ODEProblem(
        (u, p, t) -> -p[1] .* u, [1.0], (0.0, 2.0), [1.0];
        callback = ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> (i.p[1] = 3.0; nothing))
    )
    @test relerr(compare(switch, Tsit5(), 1.0e-5, switch_exact)...) < 1.0e-5
    @test relerr(compare(switch, Rosenbrock32(), 1.0e-5, switch_exact)...) < 1.0e-5
    # An affect that replaces the parameter object.
    replaced = remake(
        switch; callback = ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> (i.p = [3.0]; nothing))
    )
    @test relerr(compare(replaced, Tsit5(), 1.0e-5, switch_exact)...) < 1.0e-5

    preset = ODEProblem(
        (u, p, t) -> -p[1] .* u, [1.0], (0.0, 2.0), [1.0];
        callback = DiscreteCallback((u, t, i) -> t == 1.0, i -> (i.p[1] = 3.0; nothing)), tstops = [1.0]
    )
    @test relerr(compare(preset, Tsit5(), 1.0e-5, t -> [t < 1 ? exp(-t) : exp(-1) * exp(-3(t - 1))])...) < 1.0e-5

    # The damping changes at each zero of x: the parameters and the saltation matrix change
    # at each event.
    damped = ODEProblem(
        (u, p, t) -> [u[2], -u[1] - p[1] * u[2]], [1.0, 0.0], (0.0, 10.0), [0.1];
        callback = ContinuousCallback((u, t, i) -> u[1], i -> (i.p[1] = i.p[1] == 0.1 ? 1.0 : 0.1; nothing))
    )
    dref = solve(deepcopy(damped), Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    @test relerr(compare(damped, Tsit5(), 1.0e-5, dref)...) < 1.0e-4

    # A DAE: the algebraic equations after the event use the new parameters.
    cubic_z(x) = (q = sqrt(x^2 / 4 + 1 / 27); cbrt(x / 2 + q) + cbrt(x / 2 - q))
    cubicp!(du, u, p, t) = (du[1] = -p[1] * u[1]; du[2] = u[2]^3 + u[2] - u[1]; nothing)
    dae = ODEProblem(
        ODEFunction(cubicp!; mass_matrix = Diagonal([1.0, 0.0])), [2.0, cubic_z(2.0)], (0.0, 2.0), [1.0];
        callback = ContinuousCallback(
            (u, t, i) -> u[1] - 1.0, i -> (i.p[1] = 3.0; nothing); initializealg = BrownFullBasicInit()
        )
    )
    dae_exact = t -> (x = t < ts ? 2exp(-t) : exp(-3(t - ts)); [x, cubic_z(x)])
    fd = Rodas5P(autodiff = ADTypes.AutoFiniteDiff())
    @test relerr(compare(dae, Rodas5P(), 1.0e-5, dae_exact; adjoint_alg = fd)...) < 1.0e-4
end

@testset "GlobalAdjoint: TrajectoryError with events" begin
    # An event time is saved twice. The second entry is the error after the event,
    # e⁺ = S e⁻ with S = f⁺/f⁻ = 55 here, so the dual crosses the event first. (The first
    # entry is the linearized error before the event; the true solution has already reset
    # there, which a linearization does not see.)
    period = log(55)
    slow = ODEProblem(
        (u, p, t) -> [0.45 - u[1]], [1.0], (0.0, 6.0);
        callback = ContinuousCallback((u, t, i) -> u[1] - 0.46, nothing, i -> (i.u[1] = 1.0; nothing))
    )
    exact(t) = 0.45 + 0.55exp(-(t - fld(t, period) * period))
    sol = solve(
        slow, GlobalAdjoint(Tsit5(); gtol = 1.0, scope = TrajectoryError, samples = 1, adjoint_alg = Tsit5());
        abstol = 1.0e-5, reltol = 1.0e-5
    )
    k = findfirst(i -> sol.t[i] == sol.t[i + 1], 1:(length(sol.t) - 1))
    @test relerr(sol.global_error[k + 1], abs(exact(sol.t[k + 1]) - sol.u[k + 1][1])) < 1.0e-3
    @test sol.global_error[k + 1] > 10sol.global_error[k]
    # The entries away from the event.
    away = [j for j in 2:length(sol.t) if abs(sol.t[j] - sol.t[k]) > 0.1]
    @test all(j -> relerr(sol.global_error[j], abs(exact(sol.t[j]) - sol.u[j][1])) < 1.0e-3, away)
end

@testset "GlobalAdjoint: TrajectoryError with several calls or one-sided saves" begin
    # Doublings at t = 0.5. A callback can save before its affect, after it, or both, and
    # two callbacks at one time give four saves; each entry is the state after the calls that
    # come before it. The exact value of an entry follows from its number of doublings.
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
            ODEProblem((u, p, t) -> -u, [1.0], (0.0, 1.0)), GlobalAdjoint(Tsit5(); gtol = 1.0, scope = TrajectoryError, samples = 1, adjoint_alg = Tsit5());
            callback, tstops = [0.5], abstol = 1.0e-6, reltol = 1.0e-6
        )
        doublings = [round(Int, log2(u[1] / exp(-t))) for (t, u) in zip(sol.t, sol.u)]
        actual = [2.0^k * exp(-t) - u[1] for (k, t, u) in zip(doublings, sol.t, sol.u)]
        @test all(abs(e - abs(a)) <= 1.0e-2 * abs(a) + 1.0e-14 for (e, a) in zip(sol.global_error, actual))
    end
end

@testset "GlobalAdjoint: event edge cases" begin
    decay(u, p, t) = -p[1] .* u
    two(u, p, t) = [-p[1] * u[1], -u[2]]
    tol = 1.0e-6
    # Two DiscreteCallbacks at one time (their windows merge).
    both = ODEProblem(
        decay, [1.0], (0.0, 2.0), [1.0];
        callback = CallbackSet(
            DiscreteCallback((u, t, i) -> t == 1.0, i -> (i.u .*= 2; nothing)),
            DiscreteCallback((u, t, i) -> t == 1.0, i -> (i.p[1] = 3.0; nothing))
        ), tstops = [1.0]
    )
    @test relerr(compare(both, Tsit5(), tol, t -> [t < 1 ? exp(-t) : 2exp(-1) * exp(-3(t - 1))])...) < 1.0e-5
    # A VectorContinuousCallback with two components that have the same root.
    vcc = ODEProblem(
        two, [1.0, 1.0], (0.0, 2.0), [1.0];
        callback = VectorContinuousCallback(
            (out, u, t, i) -> (out[1] = u[1] - 0.5; out[2] = u[1] - 0.5; nothing),
            (i, idx) -> (idx == 1 ? (i.u[2] *= 2) : (i.p[1] = 3.0); nothing), 2
        )
    )
    vccref = solve(deepcopy(vcc), Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    @test relerr(compare(vcc, Tsit5(), tol, vccref)...) < 1.0e-5
    # An event in the first step (t ≈ 0.001); the reset makes it repeat.
    t1, per = log(1 / 0.999), log(2 / 0.999)
    first_step = ODEProblem(
        decay, [1.0], (0.0, 2.0), [1.0];
        callback = ContinuousCallback((u, t, i) -> u[1] - 0.999, i -> (i.u[1] = 2.0; nothing))
    )
    @test relerr(compare(first_step, Tsit5(), tol, t -> [2exp(-((t - t1) - fld(t - t1, per) * per))])...) < 1.0e-5
    # An event at the final time: the saved state is after the affect.
    at_end = ODEProblem(
        decay, [1.0], (0.0, 2.0), [1.0];
        callback = DiscreteCallback((u, t, i) -> t == 2.0, i -> (i.u .*= 2; nothing)), tstops = [2.0]
    )
    @test relerr(compare(at_end, Tsit5(), tol, t -> [2exp(-t)])...) < 1.0e-5
    # DiffEqCallbacks callbacks with state, and a SavingCallback (no change of the state).
    preset = ODEProblem(
        decay, [1.0], (0.0, 3.5), [1.0];
        callback = PresetTimeCallback([1.0, 2.0, 3.0], i -> (i.u .*= 2; nothing))
    )
    @test relerr(compare(preset, Tsit5(), tol, t -> [exp(-t) * 2.0^min(3, floor(Int, t))])...) < 1.0e-5
    periodic = ODEProblem(
        decay, [1.0], (0.0, 2.0), [1.0]; callback = PeriodicCallback(i -> (i.u .*= 1.5; nothing), 0.6)
    )
    @test relerr(compare(periodic, Tsit5(), tol, t -> [1.5^3 * exp(-t)])...) < 1.0e-5
    saving = ODEProblem(
        (u, p, t) -> [u[2], -u[1]], [1.0, 0.0], (0.0, 5.0);
        callback = SavingCallback((u, t, i) -> u[1], SavedValues(Float64, Float64))
    )
    @test relerr(compare(saving, Tsit5(), tol, t -> [cos(t), -sin(t)])...) < 1.0e-5
    # Terminating callbacks: the estimate is the error at the forward event.
    for callback in (ContinuousCallback((u, t, i) -> u[1] - 0.5, terminate!), TerminateSteadyState(1.0e-4, 1.0e-4))
        stop = ODEProblem(decay, [1.0], (0.0, 100.0), [1.0]; callback)
        @test relerr(compare(stop, Tsit5(), tol, t -> [exp(-t)])...) < 1.0e-5
    end
end

@testset "GlobalAdjoint: VJP through the events" begin
    # The default autojacvec (ReverseDiffVJP) gives Sᵀλ from AD products of the affect and the
    # condition; autojacvec = false the finite-difference saltation matrix. The estimates agree.
    estimate(prob, autojacvec; tol = 1.0e-6) = adjoint_error_estimate(
        prob, GlobalAdjoint(
            Tsit5(); samples = length(prob.u0), adjoint_alg = Tsit5(),
            sensealg = InterpolatingAdjoint(; autojacvec)
        );
        abstol = tol, reltol = tol
    )
    problems = (
        ODEProblem(
            (u, p, t) -> [1.5u[1] - u[1] * u[2], -3.0u[2] + u[1] * u[2]], [1.0, 1.0], (0.0, 10.0);
            callback = ContinuousCallback((u, t, i) -> u[1] - 3.0, i -> (i.u[1] = 1.5; nothing), nothing)
        ),
        ODEProblem(
            (u, p, t) -> [u[2], -u[1] - p[1] * u[2]], [1.0, 0.0], (0.0, 10.0), [0.1];
            callback = ContinuousCallback((u, t, i) -> u[1], i -> (i.p[1] = i.p[1] == 0.1 ? 1.0 : 0.1; nothing))
        ),
        ODEProblem(
            (u, p, t) -> [u[2], -u[1]], [0.0, 1.0], (0.0, 10.0);
            callback = VectorContinuousCallback(
                (out, u, t, i) -> (out[1] = u[1] - 0.8; out[2] = u[1] + 0.8; nothing),
                (i, idx) -> (i.u[2] = -0.9 * i.u[2]; nothing), 2
            )
        ),
        ODEProblem(
            (u, p, t) -> -u, [1.0], (0.0, 1.0);
            callback = CallbackSet(
                DiscreteCallback((u, t, i) -> t == 0.5, i -> (i.u .*= 2; nothing)),
                DiscreteCallback((u, t, i) -> t == 0.5, i -> (i.u .*= 2; nothing))
            ), tstops = [0.5]
        ),
    )
    for prob in problems
        @test relerr(estimate(prob, ReverseDiffVJP()), estimate(prob, false)) < 1.0e-5
    end

    # An affect that ForwardDiff and ReverseDiff cannot run: finite differences, with a warning.
    rigid = ODEProblem(
        (u, p, t) -> -u, [1.0], (0.0, 2.0);
        callback = ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> (i.u[1] = 2 * Float64(i.u[1]); nothing))
    )
    fallback = @test_logs (:warn, r"finite-difference saltation") match_mode = :any estimate(rigid, ReverseDiffVJP())
    # The two runs differ only in the VJP of the adjoint ODE, at the adjoint tolerance.
    @test fallback ≈ estimate(rigid, false) rtol = 1.0e-6
end

@testset "GlobalAdjoint: unsupported combinations" begin
    cb = ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> nothing)
    prob = ODEProblem((u, p, t) -> -u, [1.0], (0.0, 1.0); callback = cb)
    @test_throws ArgumentError solve(prob, GlobalAdjoint(Tsit5(); gtol = 1.0e-4, control = StepGridRefinement))
    @test_throws ArgumentError adjoint_error_estimate(prob, GlobalAdjoint(Tsit5(); sensealg = QuadratureAdjoint()))
    full = ODEProblem(
        ODEFunction((u, p, t) -> -u; mass_matrix = [1.0 0.5; 0.0 0.0]), [1.0, 1.0], (0.0, 1.0)
    )
    @test_throws ArgumentError adjoint_error_estimate(full, GlobalAdjoint(Rodas5P()))
end
