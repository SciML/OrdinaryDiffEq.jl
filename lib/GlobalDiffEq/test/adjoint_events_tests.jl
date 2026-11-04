using GlobalDiffEq, OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using SciMLSensitivity, LinearAlgebra, Test
using Random: Xoshiro
using DiffEqCallbacks: PresetTimeCallback, PeriodicCallback, SavingCallback, SavedValues,
    TerminateSteadyState
import ADTypes, SciMLBase

relerr(estimate, actual) = abs(estimate - actual) / abs(actual)

# A parameter struct without `==` (`isequal` of two copies is false).
mutable struct RateParameter
    k::Vector{Float64}
end

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
    # Global error control through the solve pipeline: the wrapper of `f` needs the dual
    # signatures for the Jacobian of Rodas5P.
    controlled = solve(smooth, GlobalAdjoint(Rodas5P(); gtol = 1.0e-4, adjoint_alg = fd))
    @test SciMLBase.successful_retcode(controlled)
    controlled_error = norm(controlled.u[end] .- [2exp(-5.0), cubic_z(2exp(-5.0))])
    @test controlled_error < 1.0e-4
    @test relerr(controlled.global_error[end], controlled_error) < 1.0e-4
    # An explicit `1.0I` is the identity.
    scaled_identity = ODEProblem(ODEFunction((u, p, t) -> -u; mass_matrix = 1.0I), [1.0, 1.0], (0.0, 1.0))
    @test adjoint_error_estimate(scaled_identity, GlobalAdjoint(Tsit5(); rng = Xoshiro(4))) ≈
        adjoint_error_estimate(ODEProblem((u, p, t) -> -u, [1.0, 1.0], (0.0, 1.0)), GlobalAdjoint(Tsit5(); rng = Xoshiro(4))) rtol = 1.0e-6
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

@testset "GlobalAdjoint: a later call at one time is linearized after the earlier ones" begin
    # At t = 1 the first callback doubles x and the second adds x² to y: the derivative of
    # the second affect is 2x at the state after the doubling (at the state before it, the
    # estimate was 0.89 of the error). Both the VJP and the finite-difference saltation.
    pair = () -> CallbackSet(
        DiscreteCallback((u, t, i) -> t == 1.0, i -> (i.u[1] *= 2; nothing)),
        DiscreteCallback((u, t, i) -> t == 1.0, i -> (i.u[2] += i.u[1]^2; nothing))
    )
    prob = ODEProblem((u, p, t) -> [-u[1], -u[2]], [1.0, 1.0], (0.0, 2.0))
    tol = 1.0e-6
    ref = solve(prob, Tsit5(); callback = pair(), tstops = [1.0], abstol = 1.0e-13, reltol = 1.0e-13)
    sol = solve(prob, Tsit5(); callback = pair(), tstops = [1.0], abstol = tol, reltol = tol)
    actual = norm(ref.u[end] .- sol.u[end])
    for sensealg in (nothing, InterpolatingAdjoint(autojacvec = false))
        ga = GlobalAdjoint(Tsit5(); samples = 2, adjoint_alg = Tsit5(), sensealg)
        estimate = adjoint_error_estimate(
            prob, ga; callback = pair(), tstops = [1.0], abstol = tol, reltol = tol
        )
        @test relerr(estimate, actual) < 1.0e-3
    end
end

@testset "GlobalAdjoint: a DiscreteCallback that a continuous event triggers" begin
    # x' = -x, y' = -y. At x = 0.5 the continuous event doubles x and sets a flag; at the
    # same step end (the event ends a step in every solve) a DiscreteCallback adds x² to y
    # and clears it. The discrete affect moves with the event time: one saltation matrix for
    # both affects (in sequence the estimate was 0.64 / 1.68 of the error).
    tied = () -> CallbackSet(
        ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> (i.u[1] *= 2; i.p[1] = 1.0; nothing)),
        DiscreteCallback((u, t, i) -> i.p[1] == 1.0, i -> (i.u[2] += i.u[1]^2; i.p[1] = 0.0; nothing))
    )
    prob = () -> ODEProblem((u, p, t) -> [-u[1], -u[2]], [1.0, 1.0], (0.0, 2.0), [0.0])
    tol = 1.0e-6
    ref = solve(prob(), Tsit5(); callback = tied(), abstol = 1.0e-13, reltol = 1.0e-13)
    sol = solve(prob(), Tsit5(); callback = tied(), abstol = tol, reltol = tol)
    actual = norm(ref.u[end] .- sol.u[end])
    for sensealg in (nothing, InterpolatingAdjoint(autojacvec = false))
        ga = GlobalAdjoint(Tsit5(); samples = 2, adjoint_alg = Tsit5(), sensealg)
        estimate = adjoint_error_estimate(prob(), ga; callback = tied(), abstol = tol, reltol = tol)
        @test relerr(estimate, actual) < 1.0e-3
    end
    # The four entries at the first event time (both callbacks save before and after their
    # affects): the linearized error against the branch of the exact solution before the
    # event, between the two affects, and after both, each continued to the event time.
    trajectory = solve(
        prob(), GlobalAdjoint(Tsit5(); gtol = 1.0, scope = TrajectoryError, samples = 2, adjoint_alg = Tsit5());
        callback = tied(), abstol = tol, reltol = tol
    )
    at_event = findall(t -> abs(t - log(2)) < 1.0e-5, trajectory.t)
    @test length(at_event) == 4
    te = trajectory.t[first(at_event)]
    branches = (
        [exp(-te), exp(-te)], [2exp(-te), exp(-te)], [2exp(-te), exp(-te)],
        [2exp(-te), exp(-te) + exp(-(te - log(2)))],
    )
    for (j, branch) in zip(at_event, branches)
        @test relerr(trajectory.global_error[j], norm(branch .- trajectory.u[j])) < 1.0e-3
    end
end

@testset "GlobalAdjoint: parameter structs without ==, a SavingCallback at an event" begin
    # The affect resets x and sets the rate k = 2. With a struct without `==`, or a vector of
    # such structs, the trial run of the affect is compared field by field (`isequal` of two
    # copies is false).
    rate(p::RateParameter) = p.k[1]
    rate(p::Vector{RateParameter}) = rate(only(p))
    rate(p) = p[1]
    set_rate!(p::RateParameter) = (p.k[1] = 2.0)
    set_rate!(p::Vector{RateParameter}) = set_rate!(only(p))
    set_rate!(p) = (p[1] = 2.0)
    reset = () -> ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> (i.u[1] = 1.0; set_rate!(i.p); nothing))
    prob = p -> ODEProblem((u, p, t) -> [-rate(p) * u[1], -u[2]], [1.0, 1.0], (0.0, 2.0), p)
    # The same directions in each estimate.
    ga = () -> GlobalAdjoint(Tsit5(); samples = 2, adjoint_alg = Tsit5(), rng = Xoshiro(5))
    tol = 1.0e-6
    vector = adjoint_error_estimate(prob([1.0]), ga(); callback = reset(), abstol = tol, reltol = tol)
    for p in (RateParameter([1.0]), [RateParameter([1.0])])
        @test adjoint_error_estimate(prob(p), ga(); callback = reset(), abstol = tol, reltol = tol) ≈
            vector rtol = 1.0e-6
    end
    # A SavingCallback with such a struct changes nothing: no event for another adjoint method.
    saving = () -> SavingCallback((u, t, i) -> copy(u), SavedValues(Float64, Vector{Float64}))
    quadrature = () -> GlobalAdjoint(Tsit5(); sensealg = QuadratureAdjoint(), rng = Xoshiro(5))
    @test adjoint_error_estimate(prob(RateParameter([1.0])), quadrature(); callback = saving()) ≈
        adjoint_error_estimate(prob(RateParameter([1.0])), quadrature()) rtol = 1.0e-6
    # The SavingCallback acts at the step end of the event too; its affect is not in the
    # saltation matrix of the event (it cannot take the dual numbers of the VJP, which then
    # used finite differences, with a warning).
    both = CallbackSet(reset(), saving())
    estimate = @test_logs min_level = Base.CoreLogging.Warn adjoint_error_estimate(
        prob([1.0]), ga(); callback = both, abstol = tol, reltol = tol
    )
    @test estimate ≈ vector rtol = 1.0e-6
end

@testset "GlobalAdjoint: affects that read the derivative, integer time span" begin
    # An affect that reads get_du is a function of u, p and t; the trial run gives f(u), which
    # differs from the forward derivative by about the local error, so the check of the trial
    # affect accepts it, and the difference is an error of the forward affect in the source
    # term (without it the estimate was 1.6 % off).
    kick = ContinuousCallback(
        (u, t, i) -> u[1], i -> (i.u[2] = -0.8 * i.u[2] + 0.1 * SciMLBase.get_du(i)[1]; nothing)
    )
    osc = ODEProblem((u, p, t) -> [u[2], -u[1]], [1.0, 0.0], (0.0, 6.0); callback = kick)
    ref = solve(osc, Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    estimate, actual = compare(osc, Tsit5(), 1.0e-7, t -> ref(t))
    @test relerr(estimate, actual) < 1.0e-4
    # An integer time span with an event.
    reset = () -> ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> (i.u[1] = 1.0; nothing))
    ga = GlobalAdjoint(Tsit5(); samples = 1, adjoint_alg = Tsit5())
    @test adjoint_error_estimate(ODEProblem((u, p, t) -> -u, [1.0], (0, 2)), ga; callback = reset()) ≈
        adjoint_error_estimate(ODEProblem((u, p, t) -> -u, [1.0], (0.0, 2.0)), ga; callback = reset()) rtol = 1.0e-6
end

@testset "GlobalAdjoint: callbacks that do not save after the affect" begin
    # Without a save after the affect, the dense output of the next step starts from the
    # state before it; the recorded solve saves the state after the event (the estimate was
    # 3 to 5 times too large).
    prob = ODEProblem((u, p, t) -> [-u[1], -u[2]], [1.0, 1.0], (0.0, 2.0))
    xexact(t) = 2.0^floor(t / log(2)) * exp(-t)
    for save_positions in ((true, true), (false, true), (true, false), (false, false))
        callback = () -> ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> (i.u[1] *= 2; nothing); save_positions)
        sol = solve(prob, Tsit5(); callback = callback(), abstol = 1.0e-7, reltol = 1.0e-7)
        estimate = adjoint_error_estimate(
            prob, GlobalAdjoint(Tsit5(); samples = 2, adjoint_alg = Tsit5());
            callback = callback(), abstol = 1.0e-7, reltol = 1.0e-7
        )
        @test relerr(estimate, norm([xexact(2.0), exp(-2.0)] .- sol.u[end])) < 1.0e-3
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

@testset "GlobalAdjoint: affect arguments, callback state, initialization, default VJP" begin
    # A VectorContinuousCallback affect gets the vector of the fired components and their
    # signs, as in the forward solve.
    walls = () -> ODEProblem(
        (u, p, t) -> [u[2], -u[1]], [0.0, 1.0], (0.0, 10.0);
        callback = VectorContinuousCallback(
            (out, u, t, i) -> (out[1] = u[1] - 0.8; out[2] = u[1] + 0.8; nothing),
            (i, ev) -> (ev[1] != 0 && (i.u[2] = -0.9 * i.u[2]); ev[2] != 0 && (i.u[2] = -0.8 * i.u[2]); nothing), 2
        )
    )
    wallref = solve(walls(), Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    @test relerr(compare(walls(), Tsit5(), 1.0e-6, wallref)...) < 1.0e-4

    # An affect with internal state: it halves the velocity only at its second call. As in
    # SciMLSensitivity, the trial runs use one copy of the callback for the whole solve, so
    # such an affect is not supported; the check finds it.
    calls = Ref(0)
    counted = ODEProblem(
        (u, p, t) -> [u[2], -u[1]], [1.0, 0.0], (0.0, 10.0);
        callback = ContinuousCallback(
            (u, t, i) -> u[1], i -> (calls[] += 1; calls[] == 2 && (i.u[2] *= 0.5); nothing)
        )
    )
    @test_throws ArgumentError adjoint_error_estimate(
        counted, GlobalAdjoint(Tsit5(); samples = 2, adjoint_alg = Tsit5()); abstol = 1.0e-6, reltol = 1.0e-6
    )
    # An affect that halves the velocity when its flag is true and flips the flag. The calls
    # that change the state all have the flag true, as the trial copy of the callback; the
    # others change nothing and are not checked. This is not support for affects with state
    # in the callback (see the counter above), but the trial runs must not change that
    # state: each one starts from its own copy of the trial callback, so the finite
    # differences of one saltation matrix all see the same flag.
    toggling = () -> (
        flag = Ref(true);
        ContinuousCallback((u, t, i) -> u[1], i -> (flag[] && (i.u[2] *= 0.5); flag[] = !flag[]; nothing))
    )
    swing = ODEProblem((u, p, t) -> [u[2], -u[1]], [1.0, 0.0], (0.0, 10.0))
    swing_ref = solve(swing, Tsit5(); callback = toggling(), abstol = 1.0e-13, reltol = 1.0e-13)
    swing_sol = solve(swing, Tsit5(); callback = toggling(), abstol = 1.0e-6, reltol = 1.0e-6)
    for sensealg in (nothing, InterpolatingAdjoint(autojacvec = false))
        estimate = adjoint_error_estimate(
            swing, GlobalAdjoint(Tsit5(); samples = 2, adjoint_alg = Tsit5(), sensealg);
            callback = toggling(), abstol = 1.0e-6, reltol = 1.0e-6
        )
        @test relerr(estimate, norm(swing_ref.u[end] .- swing_sol.u[end])) < 1.0e-4
    end

    # An affect at the initialization: the true solution starts with it too.
    initial = ODEProblem(
        (u, p, t) -> -u, [1.0], (0.0, 1.0);
        callback = PeriodicCallback(i -> (i.u .*= 2; nothing), 10.0; initial_affect = true)
    )
    @test relerr(compare(initial, Tsit5(), 1.0e-6, t -> [2exp(-t)])...) < 1.0e-5

    # A right-hand side that ReverseDiff cannot run (a Float64 cache): the default sensealg
    # uses a numerical VJP, with a warning.
    cache = zeros(2)
    cached = ODEProblem((u, p, t) -> (cache .= (-u[1], -2u[2]); copy(cache)), [1.0, 1.0], (0.0, 1.0))
    estimate, actual = @test_logs (:warn, r"numerical VJP") match_mode = :any compare(
        cached, Tsit5(), 1.0e-6, t -> [exp(-t), exp(-2t)]
    )
    @test relerr(estimate, actual) < 1.0e-4
    # The same for a DAE: the Jacobian at the final time uses finite differences, and the
    # numerical VJP an in-place right-hand side (SciMLSensitivity's finite-difference
    # Jacobian of an out-of-place DAE fails).
    dae_cached = ODEProblem(
        ODEFunction((u, p, t) -> (cache .= (-u[1], u[2] - 2u[1]); copy(cache)); mass_matrix = Diagonal([1.0, 0.0])),
        [1.0, 2.0], (0.0, 1.0)
    )
    fd = Rodas5P(autodiff = ADTypes.AutoFiniteDiff())
    estimate, actual = @test_logs (:warn, r"numerical VJP") match_mode = :any compare(
        dae_cached, fd, 1.0e-6, t -> [exp(-t), 2exp(-t)]; adjoint_alg = fd
    )
    @test relerr(estimate, actual) < 1.0e-3
end

@testset "GlobalAdjoint: unsupported combinations" begin
    cb = ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> nothing)
    prob = ODEProblem((u, p, t) -> -u, [1.0], (0.0, 1.0); callback = cb)
    @test_throws ArgumentError solve(prob, GlobalAdjoint(Tsit5(); gtol = 1.0e-4, control = StepGridRefinement))
    # Another adjoint method: an event that changes the state, or a DAE, fails (the event at
    # its call, the DAE before the solve); a call that changes nothing (a SavingCallback) is
    # no event.
    quadrature = GlobalAdjoint(Tsit5(); sensealg = QuadratureAdjoint(), samples = 1)
    reset = ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> (i.u[1] = 1.0; nothing))
    decay = ODEProblem((u, p, t) -> -u, [1.0], (0.0, 1.0))
    @test_throws ArgumentError adjoint_error_estimate(decay, quadrature; callback = reset)
    @test_throws ArgumentError adjoint_error_estimate(
        ODEProblem(ODEFunction((u, p, t) -> [-u[1], u[2] - u[1]]; mass_matrix = Diagonal([1.0, 0.0])), [1.0, 1.0], (0.0, 1.0)),
        GlobalAdjoint(Rodas5P(); sensealg = QuadratureAdjoint())
    )
    saving = SavingCallback((u, t, i) -> copy(u), SavedValues(Float64, Vector{Float64}))
    @test adjoint_error_estimate(decay, quadrature; callback = saving) ≈
        adjoint_error_estimate(decay, quadrature) rtol = 1.0e-6
    @test_throws ArgumentError adjoint_error_estimate(
        prob, GlobalAdjoint(Tsit5(); sensealg = InterpolatingAdjoint(checkpointing = true))
    )
    full = ODEProblem(
        ODEFunction((u, p, t) -> -u; mass_matrix = [1.0 0.5; 0.0 0.0]), [1.0, 1.0], (0.0, 1.0)
    )
    @test_throws ArgumentError adjoint_error_estimate(full, GlobalAdjoint(Rodas5P()))
    # A multiple of the identity is not the identity: the error identity assumes M = I.
    scaled = ODEProblem(ODEFunction((u, p, t) -> -u; mass_matrix = 2.0I), [1.0, 1.0], (0.0, 1.0))
    @test_throws ArgumentError adjoint_error_estimate(scaled, GlobalAdjoint(Rodas5P()))
    @test_throws ArgumentError solve(scaled, GlobalAdjoint(Rodas5P(); gtol = 1.0e-4))
end
