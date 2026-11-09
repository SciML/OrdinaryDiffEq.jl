using GlobalDiffEq, OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve
using SciMLSensitivity, LinearAlgebra, Test
using Random: Xoshiro
using DiffEqCallbacks: PresetTimeCallback, PeriodicCallback, SavingCallback, SavedValues,
    TerminateSteadyState
import ADTypes, SciMLBase, SciMLLogging

relerr(estimate, actual) = abs(estimate - actual) / abs(actual)

# Return the estimate and the actual error of a plain solve at the final time. One
# projection per state gives the 2-norm of the linearized error, with no sampling error.
# The plain solve uses a copy because an affect can change the parameters in place.
function compare(prob, alg, tol, exact; adjoint_alg = Tsit5())
    ga = GlobalAdjoint(alg; samples = length(prob.u0), adjoint_alg)
    estimate = adjoint_error_estimate(prob, ga; abstol = tol, reltol = tol)
    sol = solve(deepcopy(prob), alg; abstol = tol, reltol = tol)
    return estimate, norm(exact(sol.t[end]) .- sol.u[end])
end

# A parameter struct without `==` (`isequal` of two copies is false).
mutable struct RateParameter
    k::Vector{Float64}
end
rate(p::RateParameter) = p.k[1]
rate(p::Vector{RateParameter}) = rate(only(p))
rate(p) = p[1]
set_rate!(p::RateParameter) = (p.k[1] = 2.0)
set_rate!(p::Vector{RateParameter}) = set_rate!(only(p))
set_rate!(p) = (p[1] = 2.0)

# The tests use named right-hand sides, conditions and affects. Each new function type
# compiles the solver and the adjoint stack again, so the tests share these functions.
decay(u, p, t) = -p[1] .* u
decay_pair(u, p, t) = [-u[1], -u[2]]
rate_decay(u, p, t) = [-rate(p) * u[1], -u[2]]
oscillator(u, p, t) = [u[2], -u[1]]
damped(u, p, t) = [u[2], -u[1] - p[1] * u[2]]
lotka_volterra(u, p, t) = [1.5u[1] - u[1] * u[2], -3.0u[2] + u[1] * u[2]]
relax(u, p, t) = [0.45 - u[1]]
# A DAE with the algebraic equation z³ + z = x, and its solution z(x).
cubic!(du, u, p, t) = (du[1] = -p[1] * u[1]; du[2] = u[2]^3 + u[2] - u[1]; nothing)
cubic_z(x) = (q = sqrt(x^2 / 4 + 1 / 27); cbrt(x / 2 + q) + cbrt(x / 2 - q))
cubic_dae(u0, tspan; kwargs...) =
    ODEProblem(ODEFunction(cubic!; mass_matrix = Diagonal([1.0, 0.0])), u0, tspan, [1.0]; kwargs...)
# ReverseDiff cannot run these right-hand sides, because they write to a Float64 cache.
const CACHE = zeros(2)
cached_decay(u, p, t) = (CACHE .= (-u[1], -2u[2]); copy(CACHE))
cached_dae(u, p, t) = (CACHE .= (-u[1], u[2] - 2u[1]); copy(CACHE))

below_half(u, t, i) = u[1] - 0.5
double!(i) = (i.u .*= 2; nothing)
fd = Rodas5P(autodiff = ADTypes.AutoFiniteDiff())

# Problems that more than one test set uses. Each call gives new parameters, because some
# affects change the parameters in place.
reset_callback() = ContinuousCallback(below_half, i -> (i.u[1] = 1.0; nothing))
saving_callback() = SavingCallback((u, t, i) -> copy(u), SavedValues(Float64, Vector{Float64}))
lotka_volterra_problem() = ODEProblem(
    lotka_volterra, [1.0, 1.0], (0.0, 10.0);
    callback = ContinuousCallback((u, t, i) -> u[1] - 3.0, i -> (i.u[1] = 1.5; nothing), nothing)
)
# The damping switches at each zero of x.
damped_problem() = ODEProblem(
    damped, [1.0, 0.0], (0.0, 10.0), [0.1];
    callback = ContinuousCallback((u, t, i) -> u[1], i -> (i.p[1] = i.p[1] == 0.1 ? 1.0 : 0.1; nothing))
)
walls_problem() = ODEProblem(
    oscillator, [0.0, 1.0], (0.0, 10.0);
    callback = VectorContinuousCallback(
        (out, u, t, i) -> (out[1] = u[1] - 0.8; out[2] = u[1] + 0.8; nothing),
        (i, idx) -> (i.u[2] = -0.9 * i.u[2]; nothing), 2
    )
)
# x relaxes to 0.45 and resets to 1 at 0.46, with the period log(55).
relax_problem(t_end) = ODEProblem(
    relax, [1.0], (0.0, t_end);
    callback = ContinuousCallback((u, t, i) -> u[1] - 0.46, nothing, i -> (i.u[1] = 1.0; nothing))
)
relax_exact(t) = (period = log(55); 0.45 + 0.55exp(-(t - fld(t, period) * period)))
# A GlobalAdjoint that returns the error of each saved state.
trajectory_alg(samples) = GlobalAdjoint(Tsit5(); gtol = 1.0, scope = TrajectoryError, samples, adjoint_alg = Tsit5())
# Two DiscreteCallbacks that double u at t = 0.5.
doubling_pair() = CallbackSet(
    DiscreteCallback((u, t, i) -> t == 0.5, double!), DiscreteCallback((u, t, i) -> t == 0.5, double!)
)

@testset "GlobalAdjoint: jumps of the dense output" begin
    # The dense output of Rosenbrock32 does not end at the state of the step, so the error
    # identity needs the knot terms and the endpoint shift.
    prob = ODEProblem(decay, [1.0], (0.0, 3.0), [1.0])
    @test relerr(compare(prob, Rosenbrock32(), 1.0e-5, t -> [exp(-t)])...) < 1.0e-5
end

@testset "GlobalAdjoint: events" begin
    lv = lotka_volterra_problem()
    lvref = solve(lv, Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    # The estimate is a linearization, so its relative error falls with the tolerance.
    @test relerr(compare(lv, Tsit5(), 1.0e-3, lvref)...) < 1.0e-2
    @test relerr(compare(lv, Tsit5(), 1.0e-5, lvref)...) < 1.0e-3

    @test relerr(compare(relax_problem(10.0), Tsit5(), 1.0e-3, t -> [relax_exact(t)])...) < 1.0e-3

    preset = ODEProblem(
        decay, [1.0], (0.0, 3.5), [1.0];
        callback = DiscreteCallback((u, t, i) -> t in (1.0, 2.0, 3.0), double!), tstops = [1.0, 2.0, 3.0]
    )
    @test relerr(compare(preset, Tsit5(), 1.0e-4, t -> [exp(-t) * 2.0^min(3, floor(Int, t))])...) < 1.0e-4

    walls = walls_problem()
    wallref = solve(walls, Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    @test relerr(compare(walls, Tsit5(), 1.0e-5, wallref)...) < 1.0e-3

    # Global error control with events.
    controlled = solve(lv, GlobalAdjoint(Tsit5(); gtol = 1.0e-5, samples = 2, adjoint_alg = Tsit5()))
    @test norm(lvref(10.0) .- controlled.u[end]) < 2.0e-5
end

@testset "GlobalAdjoint: DAEs" begin
    smooth = cubic_dae([2.0, 1.0], (0.0, 5.0))
    estimate, actual = compare(smooth, Rodas5P(), 1.0e-4, t -> (x = 2exp(-t); [x, cubic_z(x)]); adjoint_alg = fd)
    @test relerr(estimate, actual) < 1.0e-4
    per = log(4.0)
    reset = cubic_dae(
        [2.0, 1.0], (0.0, 5.0);
        callback = ContinuousCallback(
            (u, t, i) -> u[1] - 0.5, nothing, i -> (i.u[1] = 2.0; nothing);
            initializealg = BrownFullBasicInit()
        )
    )
    exact = t -> (x = 2exp(-(t - fld(t, per) * per)); [x, cubic_z(x)])
    @test relerr(compare(reset, Rodas5P(), 1.0e-4, exact; adjoint_alg = fd)...) < 1.0e-3
    # Global error control through `solve`. The wrapper of `f` needs the dual signatures
    # for the Jacobian of Rodas5P.
    controlled = solve(smooth, GlobalAdjoint(Rodas5P(); gtol = 1.0e-4, adjoint_alg = fd))
    @test SciMLBase.successful_retcode(controlled)
    controlled_error = norm(controlled.u[end] .- [2exp(-5.0), cubic_z(2exp(-5.0))])
    @test controlled_error < 1.0e-4
    @test relerr(controlled.global_error[end], controlled_error) < 1.0e-4
    # A mass matrix of `1.0I` counts as the identity.
    scaled_identity = ODEProblem(ODEFunction(decay; mass_matrix = 1.0I), [1.0, 1.0], (0.0, 1.0), [1.0])
    plain = ODEProblem(decay, [1.0, 1.0], (0.0, 1.0), [1.0])
    @test adjoint_error_estimate(scaled_identity, GlobalAdjoint(Tsit5(); rng = Xoshiro(4))) ≈
        adjoint_error_estimate(plain, GlobalAdjoint(Tsit5(); rng = Xoshiro(4))) rtol = 1.0e-6
end

@testset "GlobalAdjoint: events that change the parameters" begin
    # The rate changes from 1 to 3 when u reaches 1/2. The dual uses the parameters of each
    # interval. The saltation matrix uses the vector field after the event (S = φ⁺/φ⁻ = 3).
    ts = log(2.0)
    switch_exact = t -> [t < ts ? exp(-t) : 0.5exp(-3(t - ts))]
    switch = ODEProblem(
        decay, [1.0], (0.0, 2.0), [1.0]; callback = ContinuousCallback(below_half, i -> (i.p[1] = 3.0; nothing))
    )
    @test relerr(compare(switch, Tsit5(), 1.0e-5, switch_exact)...) < 1.0e-5
    @test relerr(compare(switch, Rosenbrock32(), 1.0e-5, switch_exact)...) < 1.0e-5
    # An affect that replaces the parameter object.
    replaced = remake(switch; callback = ContinuousCallback(below_half, i -> (i.p = [3.0]; nothing)))
    @test relerr(compare(replaced, Tsit5(), 1.0e-5, switch_exact)...) < 1.0e-5

    preset = ODEProblem(
        decay, [1.0], (0.0, 2.0), [1.0];
        callback = DiscreteCallback((u, t, i) -> t == 1.0, i -> (i.p[1] = 3.0; nothing)), tstops = [1.0]
    )
    @test relerr(compare(preset, Tsit5(), 1.0e-5, t -> [t < 1 ? exp(-t) : exp(-1) * exp(-3(t - 1))])...) < 1.0e-5

    # The parameters and the saltation matrix change at each event.
    damped = damped_problem()
    dref = solve(deepcopy(damped), Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    @test relerr(compare(damped, Tsit5(), 1.0e-5, dref)...) < 1.0e-4

    # A DAE whose rate changes at the event.
    dae = cubic_dae(
        [2.0, cubic_z(2.0)], (0.0, 2.0);
        callback = ContinuousCallback(
            (u, t, i) -> u[1] - 1.0, i -> (i.p[1] = 3.0; nothing); initializealg = BrownFullBasicInit()
        )
    )
    dae_exact = t -> (x = t < ts ? 2exp(-t) : exp(-3(t - ts)); [x, cubic_z(x)])
    @test relerr(compare(dae, Rodas5P(), 1.0e-5, dae_exact; adjoint_alg = fd)...) < 1.0e-4
end

@testset "GlobalAdjoint: TrajectoryError with events" begin
    # The solution saves each event time twice. The second entry is the error after the
    # event, e⁺ = S e⁻ with S = φ⁺/φ⁻ = 55, so its dual must cross the event. The first
    # entry is the linearized error before the event, which cannot see the true reset.
    sol = solve(relax_problem(6.0), trajectory_alg(1); abstol = 1.0e-5, reltol = 1.0e-5)
    k = findfirst(i -> sol.t[i] == sol.t[i + 1], 1:(length(sol.t) - 1))
    @test relerr(sol.global_error[k + 1], abs(relax_exact(sol.t[k + 1]) - sol.u[k + 1][1])) < 1.0e-3
    @test sol.global_error[k + 1] > 10sol.global_error[k]
    # The entries away from the event match the actual error.
    away = [j for j in 2:length(sol.t) if abs(sol.t[j] - sol.t[k]) > 0.1]
    @test all(j -> relerr(sol.global_error[j], abs(relax_exact(sol.t[j]) - sol.u[j][1])) < 1.0e-3, away)
end

@testset "GlobalAdjoint: TrajectoryError with several calls or one-sided saves" begin
    # Each case doubles u at t = 0.5. A callback can save before its affect, after it, or
    # both; two callbacks at one time give four saves. Each entry is the state after the
    # calls before it, so its exact value follows from its number of doublings.
    cases = (
        DiscreteCallback((u, t, i) -> t == 0.5, double!; save_positions = (true, false)),
        DiscreteCallback((u, t, i) -> t == 0.5, double!; save_positions = (false, true)),
        doubling_pair(),
    )
    for callback in cases
        sol = solve(
            ODEProblem(decay, [1.0], (0.0, 1.0), [1.0]), trajectory_alg(1);
            callback, tstops = [0.5], abstol = 1.0e-6, reltol = 1.0e-6
        )
        doublings = [round(Int, log2(u[1] / exp(-t))) for (t, u) in zip(sol.t, sol.u)]
        actual = [2.0^k * exp(-t) - u[1] for (k, t, u) in zip(doublings, sol.t, sol.u)]
        @test all(abs(e - abs(a)) <= 1.0e-2 * abs(a) + 1.0e-14 for (e, a) in zip(sol.global_error, actual))
    end
end

@testset "GlobalAdjoint: a later call at one time is linearized after the earlier ones" begin
    # At t = 1 the first callback doubles x and the second adds x² to y. The derivative of
    # the second affect must be 2x at the state after the doubling, not before it. The test
    # checks both the VJP and the finite-difference saltation matrix.
    pair = () -> CallbackSet(
        DiscreteCallback((u, t, i) -> t == 1.0, i -> (i.u[1] *= 2; nothing)),
        DiscreteCallback((u, t, i) -> t == 1.0, i -> (i.u[2] += i.u[1]^2; nothing))
    )
    prob = ODEProblem(decay_pair, [1.0, 1.0], (0.0, 2.0))
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
    # x' = -x, y' = -y. At x = 0.5 the continuous event doubles x and sets a flag. At the
    # same step end, a DiscreteCallback adds x² to y and clears the flag. The discrete
    # affect moves with the event time, so both affects need one shared saltation matrix.
    tied = () -> CallbackSet(
        ContinuousCallback(below_half, i -> (i.u[1] *= 2; i.p[1] = 1.0; nothing)),
        DiscreteCallback((u, t, i) -> i.p[1] == 1.0, i -> (i.u[2] += i.u[1]^2; i.p[1] = 0.0; nothing))
    )
    prob = () -> ODEProblem(decay_pair, [1.0, 1.0], (0.0, 2.0), [0.0])
    tol = 1.0e-6
    ref = solve(prob(), Tsit5(); callback = tied(), abstol = 1.0e-13, reltol = 1.0e-13)
    sol = solve(prob(), Tsit5(); callback = tied(), abstol = tol, reltol = tol)
    actual = norm(ref.u[end] .- sol.u[end])
    for sensealg in (nothing, InterpolatingAdjoint(autojacvec = false))
        ga = GlobalAdjoint(Tsit5(); samples = 2, adjoint_alg = Tsit5(), sensealg)
        estimate = adjoint_error_estimate(prob(), ga; callback = tied(), abstol = tol, reltol = tol)
        @test relerr(estimate, actual) < 1.0e-3
    end
    # Both callbacks save before and after their affects, so the first event time has four
    # entries. Each entry matches the exact branch before the event, between the two
    # affects, or after both, continued to the event time.
    trajectory = solve(prob(), trajectory_alg(2); callback = tied(), abstol = tol, reltol = tol)
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
    # The affect resets x and sets the rate k = 2. The check of the trial affect compares
    # the parameters field by field, because `isequal` of two copies of a struct without
    # `==` is false. The test uses such a struct and a vector of them.
    reset = () -> ContinuousCallback(below_half, i -> (i.u[1] = 1.0; set_rate!(i.p); nothing))
    prob = p -> ODEProblem(rate_decay, [1.0, 1.0], (0.0, 2.0), p)
    # A fixed rng gives the same directions in each estimate.
    ga = () -> GlobalAdjoint(Tsit5(); samples = 2, adjoint_alg = Tsit5(), rng = Xoshiro(5))
    tol = 1.0e-6
    vector = adjoint_error_estimate(prob([1.0]), ga(); callback = reset(), abstol = tol, reltol = tol)
    for p in (RateParameter([1.0]), [RateParameter([1.0])])
        @test adjoint_error_estimate(prob(p), ga(); callback = reset(), abstol = tol, reltol = tol) ≈
            vector rtol = 1.0e-6
    end
    # With another adjoint method, a SavingCallback with such a struct is not an event.
    quadrature = () -> GlobalAdjoint(Tsit5(); sensealg = QuadratureAdjoint(), rng = Xoshiro(5))
    @test adjoint_error_estimate(prob(RateParameter([1.0])), quadrature(); callback = saving_callback()) ≈
        adjoint_error_estimate(prob(RateParameter([1.0])), quadrature()) rtol = 1.0e-6
    # The SavingCallback also acts at the step end of the event. Its affect must stay out
    # of the saltation matrix: it cannot take the dual numbers of the VJP.
    both = CallbackSet(reset(), saving_callback())
    estimate = @test_logs min_level = Base.CoreLogging.Warn adjoint_error_estimate(
        prob([1.0]), ga(); callback = both, abstol = tol, reltol = tol
    )
    @test estimate ≈ vector rtol = 1.0e-6
end

@testset "GlobalAdjoint: affects that read the derivative, integer time span" begin
    # In the trial run, get_du gives f(u), which differs from the forward derivative by
    # about the local error. The check of the trial affect accepts this, and the source
    # term counts the difference as an error of the forward affect.
    kick = ContinuousCallback(
        (u, t, i) -> u[1], i -> (i.u[2] = -0.8 * i.u[2] + 0.1 * SciMLBase.get_du(i)[1]; nothing)
    )
    osc = ODEProblem(oscillator, [1.0, 0.0], (0.0, 6.0); callback = kick)
    ref = solve(osc, Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    estimate, actual = compare(osc, Tsit5(), 1.0e-7, t -> ref(t))
    @test relerr(estimate, actual) < 1.0e-4
    # An integer time span with an event gives the same estimate as a float one.
    ga = GlobalAdjoint(Tsit5(); samples = 1, adjoint_alg = Tsit5())
    integer_span = ODEProblem(decay, [1.0], (0, 2), [1.0])
    float_span = ODEProblem(decay, [1.0], (0.0, 2.0), [1.0])
    @test adjoint_error_estimate(integer_span, ga; callback = reset_callback()) ≈
        adjoint_error_estimate(float_span, ga; callback = reset_callback()) rtol = 1.0e-6
end

@testset "GlobalAdjoint: callbacks that do not save after the affect" begin
    # Without a save after the affect, the dense output of the next step starts from the
    # state before the affect. The recorded solve must save the state after the event.
    prob = ODEProblem(decay_pair, [1.0, 1.0], (0.0, 2.0))
    xexact(t) = 2.0^floor(t / log(2)) * exp(-t)
    for save_positions in ((true, true), (false, true), (true, false), (false, false))
        callback = () -> ContinuousCallback(below_half, i -> (i.u[1] *= 2; nothing); save_positions)
        sol = solve(prob, Tsit5(); callback = callback(), abstol = 1.0e-7, reltol = 1.0e-7)
        estimate = adjoint_error_estimate(
            prob, GlobalAdjoint(Tsit5(); samples = 2, adjoint_alg = Tsit5());
            callback = callback(), abstol = 1.0e-7, reltol = 1.0e-7
        )
        @test relerr(estimate, norm([xexact(2.0), exp(-2.0)] .- sol.u[end])) < 1.0e-3
    end
end

@testset "GlobalAdjoint: event edge cases" begin
    tol = 1.0e-6
    # Two DiscreteCallbacks at one time (their windows merge).
    both = ODEProblem(
        decay, [1.0], (0.0, 2.0), [1.0];
        callback = CallbackSet(
            DiscreteCallback((u, t, i) -> t == 1.0, double!),
            DiscreteCallback((u, t, i) -> t == 1.0, i -> (i.p[1] = 3.0; nothing))
        ), tstops = [1.0]
    )
    @test relerr(compare(both, Tsit5(), tol, t -> [t < 1 ? exp(-t) : 2exp(-1) * exp(-3(t - 1))])...) < 1.0e-5
    # A VectorContinuousCallback with two components that have the same root.
    vcc = ODEProblem(
        rate_decay, [1.0, 1.0], (0.0, 2.0), [1.0];
        callback = VectorContinuousCallback(
            (out, u, t, i) -> (out[1] = u[1] - 0.5; out[2] = u[1] - 0.5; nothing),
            (i, idx) -> (idx == 1 ? (i.u[2] *= 2) : (i.p[1] = 3.0); nothing), 2
        )
    )
    vccref = solve(deepcopy(vcc), Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    @test relerr(compare(vcc, Tsit5(), tol, vccref)...) < 1.0e-5
    # An event in the first step (t ≈ 0.001). The reset makes the event repeat.
    t1, per = log(1 / 0.999), log(2 / 0.999)
    first_step = ODEProblem(
        decay, [1.0], (0.0, 2.0), [1.0];
        callback = ContinuousCallback((u, t, i) -> u[1] - 0.999, i -> (i.u[1] = 2.0; nothing))
    )
    @test relerr(compare(first_step, Tsit5(), tol, t -> [2exp(-((t - t1) - fld(t - t1, per) * per))])...) < 1.0e-5
    # An event at the final time: the saved state is after the affect.
    at_end = ODEProblem(
        decay, [1.0], (0.0, 2.0), [1.0];
        callback = DiscreteCallback((u, t, i) -> t == 2.0, double!), tstops = [2.0]
    )
    @test relerr(compare(at_end, Tsit5(), tol, t -> [2exp(-t)])...) < 1.0e-5
    # DiffEqCallbacks callbacks with internal state, and a SavingCallback that does not
    # change the state.
    preset = ODEProblem(decay, [1.0], (0.0, 3.5), [1.0]; callback = PresetTimeCallback([1.0, 2.0, 3.0], double!))
    @test relerr(compare(preset, Tsit5(), tol, t -> [exp(-t) * 2.0^min(3, floor(Int, t))])...) < 1.0e-5
    periodic = ODEProblem(
        decay, [1.0], (0.0, 2.0), [1.0]; callback = PeriodicCallback(i -> (i.u .*= 1.5; nothing), 0.6)
    )
    @test relerr(compare(periodic, Tsit5(), tol, t -> [1.5^3 * exp(-t)])...) < 1.0e-5
    saving = ODEProblem(
        oscillator, [1.0, 0.0], (0.0, 5.0);
        callback = SavingCallback((u, t, i) -> u[1], SavedValues(Float64, Float64))
    )
    @test relerr(compare(saving, Tsit5(), tol, t -> [cos(t), -sin(t)])...) < 1.0e-5
    # Terminating callbacks: the estimate is the error at the forward event.
    for callback in (ContinuousCallback(below_half, terminate!), TerminateSteadyState(1.0e-4, 1.0e-4))
        stop = ODEProblem(decay, [1.0], (0.0, 100.0), [1.0]; callback)
        @test relerr(compare(stop, Tsit5(), tol, t -> [exp(-t)])...) < 1.0e-5
    end
end

@testset "GlobalAdjoint: VJP through the events" begin
    # The default autojacvec (ReverseDiffVJP) computes Sᵀλ from AD products of the affect
    # and the condition. autojacvec = false uses the finite-difference saltation matrix.
    # The two estimates agree.
    estimate(prob, autojacvec; tol = 1.0e-6) = adjoint_error_estimate(
        prob, GlobalAdjoint(
            Tsit5(); samples = length(prob.u0), adjoint_alg = Tsit5(),
            sensealg = InterpolatingAdjoint(; autojacvec)
        );
        abstol = tol, reltol = tol
    )
    problems = (
        lotka_volterra_problem(), damped_problem(), walls_problem(),
        ODEProblem(decay, [1.0], (0.0, 1.0), [1.0]; callback = doubling_pair(), tstops = [0.5]),
    )
    for prob in problems
        @test relerr(estimate(prob, ReverseDiffVJP()), estimate(prob, false)) < 1.0e-5
    end

    # ForwardDiff and ReverseDiff cannot run this affect, so the saltation matrix uses
    # finite differences and a warning says so.
    rigid = ODEProblem(
        decay, [1.0], (0.0, 2.0), [1.0];
        callback = ContinuousCallback(below_half, i -> (i.u[1] = 2 * Float64(i.u[1]); nothing))
    )
    fallback = @test_logs (:warn, r"finite-difference saltation") match_mode = :any estimate(rigid, ReverseDiffVJP())
    # The two runs differ only in the VJP of the adjoint ODE, so they agree to the adjoint
    # tolerance.
    @test fallback ≈ estimate(rigid, false) rtol = 1.0e-6
end

@testset "GlobalAdjoint: VJP through the events of a DAE" begin
    # The VJP form goes through the algebraic projection after each affect with the implicit
    # function theorem. It agrees with the finite-difference saltation matrix. Without a
    # warning, no event fell back to finite differences.
    @test GlobalDiffEq._algebraic_solver() !== nothing
    estimate(prob, autojacvec) = adjoint_error_estimate(
        prob, GlobalAdjoint(
            Rodas5P(); samples = 1, adjoint_alg = fd, rng = Xoshiro(6),
            sensealg = InterpolatingAdjoint(; autojacvec)
        );
        abstol = 1.0e-6, reltol = 1.0e-6
    )
    problems = (
        # A reset of x, with the reinitialization of z.
        cubic_dae(
            [2.0, 1.0], (0.0, 3.0);
            callback = ContinuousCallback(
                below_half, nothing, i -> (i.u[1] = 2.0; nothing); initializealg = BrownFullBasicInit()
            )
        ),
        # The rate changes at the event.
        cubic_dae(
            [2.0, cubic_z(2.0)], (0.0, 2.0);
            callback = ContinuousCallback(
                (u, t, i) -> u[1] - 1.0, i -> (i.p[1] = 3.0; nothing); initializealg = BrownFullBasicInit()
            )
        ),
        # Two calls at one time: the second sees the state after the first.
        cubic_dae(
            [2.0, cubic_z(2.0)], (0.0, 2.0);
            callback = CallbackSet(
                DiscreteCallback((u, t, i) -> t == 1.0, double!; initializealg = BrownFullBasicInit()),
                DiscreteCallback(
                    (u, t, i) -> t == 1.0, i -> (i.p[1] = 3.0; nothing); initializealg = BrownFullBasicInit()
                )
            ), tstops = [1.0]
        ),
        # The affect reads the algebraic variable, so A contains dz/dx.
        cubic_dae(
            [2.0, cubic_z(2.0)], (0.0, 2.0);
            callback = ContinuousCallback(
                below_half, i -> (i.u[1] = 1.0 + i.u[2]; nothing); initializealg = BrownFullBasicInit()
            )
        ),
    )
    for prob in problems
        vjp = @test_logs min_level = Base.CoreLogging.Warn estimate(prob, ReverseDiffVJP())
        @test relerr(vjp, estimate(prob, false)) < 1.0e-5
    end
end

@testset "GlobalAdjoint: affect arguments, callback state, initialization, default VJP" begin
    # The affect of a VectorContinuousCallback gets the fired components with their signs,
    # as in the forward solve.
    walls = () -> ODEProblem(
        oscillator, [0.0, 1.0], (0.0, 10.0);
        callback = VectorContinuousCallback(
            (out, u, t, i) -> (out[1] = u[1] - 0.8; out[2] = u[1] + 0.8; nothing),
            (i, ev) -> (ev[1] != 0 && (i.u[2] = -0.9 * i.u[2]); ev[2] != 0 && (i.u[2] = -0.8 * i.u[2]); nothing), 2
        )
    )
    wallref = solve(walls(), Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    @test relerr(compare(walls(), Tsit5(), 1.0e-6, wallref)...) < 1.0e-4

    # This affect has internal state: it halves the velocity only at its second call. As in
    # SciMLSensitivity, the trial runs use one copy of the callback for the whole solve, so
    # the check rejects such an affect.
    calls = Ref(0)
    counted = ODEProblem(
        oscillator, [1.0, 0.0], (0.0, 10.0);
        callback = ContinuousCallback(
            (u, t, i) -> u[1], i -> (calls[] += 1; calls[] == 2 && (i.u[2] *= 0.5); nothing)
        )
    )
    @test_throws ArgumentError adjoint_error_estimate(
        counted, GlobalAdjoint(Tsit5(); samples = 2, adjoint_alg = Tsit5()); abstol = 1.0e-6, reltol = 1.0e-6
    )
    # This affect halves the velocity when its flag is true, then flips the flag. The calls
    # that change the state all see the flag true, as the trial copy does, so the check
    # passes. Each trial run starts from its own copy, so all finite differences of one
    # saltation matrix see the same flag.
    toggling = () -> (
        flag = Ref(true);
        ContinuousCallback((u, t, i) -> u[1], i -> (flag[] && (i.u[2] *= 0.5); flag[] = !flag[]; nothing))
    )
    swing = ODEProblem(oscillator, [1.0, 0.0], (0.0, 10.0))
    swing_ref = solve(swing, Tsit5(); callback = toggling(), abstol = 1.0e-13, reltol = 1.0e-13)
    swing_sol = solve(swing, Tsit5(); callback = toggling(), abstol = 1.0e-6, reltol = 1.0e-6)
    for sensealg in (nothing, InterpolatingAdjoint(autojacvec = false))
        estimate = adjoint_error_estimate(
            swing, GlobalAdjoint(Tsit5(); samples = 2, adjoint_alg = Tsit5(), sensealg);
            callback = toggling(), abstol = 1.0e-6, reltol = 1.0e-6
        )
        @test relerr(estimate, norm(swing_ref.u[end] .- swing_sol.u[end])) < 1.0e-4
    end

    # An affect at initialization: the true solution starts with the same affect.
    initial = ODEProblem(
        decay, [1.0], (0.0, 1.0), [1.0]; callback = PeriodicCallback(double!, 10.0; initial_affect = true)
    )
    @test relerr(compare(initial, Tsit5(), 1.0e-6, t -> [2exp(-t)])...) < 1.0e-5

    # ReverseDiff cannot run this right-hand side, so the default sensealg uses a
    # numerical VJP and warns.
    cached = ODEProblem(cached_decay, [1.0, 1.0], (0.0, 1.0))
    estimate, actual = @test_logs (:warn, r"numerical VJP") match_mode = :any compare(
        cached, Tsit5(), 1.0e-6, t -> [exp(-t), exp(-2t)]
    )
    @test relerr(estimate, actual) < 1.0e-4
    # The `global_error_ad_fallback` toggle of `DEVerbosity` silences the warning.
    @test_logs min_level = Base.CoreLogging.Warn adjoint_error_estimate(
        cached, GlobalAdjoint(Tsit5()); abstol = 1.0e-6, reltol = 1.0e-6,
        verbose = DEVerbosity(global_error_ad_fallback = SciMLLogging.Silent())
    )
    # The same fallback for a DAE. The Jacobian at the final time uses finite differences,
    # and the numerical VJP uses an in-place right-hand side. The finite-difference Jacobian
    # of SciMLSensitivity fails for an out-of-place DAE.
    dae_cached = ODEProblem(ODEFunction(cached_dae; mass_matrix = Diagonal([1.0, 0.0])), [1.0, 2.0], (0.0, 1.0))
    estimate, actual = @test_logs (:warn, r"numerical VJP") match_mode = :any compare(
        dae_cached, fd, 1.0e-6, t -> [exp(-t), 2exp(-t)]; adjoint_alg = fd
    )
    @test relerr(estimate, actual) < 1.0e-3
end

@testset "GlobalAdjoint: the message at a grazing event" begin
    # At a grazing event the saltation matrix has no time-shift term, and a warning says so.
    # The `global_error_grazing_event` toggle of `DEVerbosity` silences it.
    @test_logs (:warn, r"t = 1.5 is grazing") GlobalDiffEq._grazing_message(1.5, DEVerbosity())
    @test_logs GlobalDiffEq._grazing_message(
        1.5, DEVerbosity(global_error_grazing_event = SciMLLogging.Silent())
    )
end

@testset "GlobalAdjoint: unsupported combinations" begin
    prob = ODEProblem(decay, [1.0], (0.0, 1.0), [1.0]; callback = ContinuousCallback(below_half, i -> nothing))
    @test_throws ArgumentError solve(prob, GlobalAdjoint(Tsit5(); gtol = 1.0e-4, control = StepGridRefinement))
    # With another adjoint method, an event that changes the state fails at its call, and a
    # DAE fails before the solve. A SavingCallback changes nothing, so it is not an event.
    quadrature = GlobalAdjoint(Tsit5(); sensealg = QuadratureAdjoint(), samples = 1)
    plain = ODEProblem(decay, [1.0], (0.0, 1.0), [1.0])
    @test_throws ArgumentError adjoint_error_estimate(plain, quadrature; callback = reset_callback())
    @test_throws ArgumentError adjoint_error_estimate(
        cubic_dae([2.0, 1.0], (0.0, 1.0)), GlobalAdjoint(Rodas5P(); sensealg = QuadratureAdjoint())
    )
    @test adjoint_error_estimate(plain, quadrature; callback = saving_callback()) ≈
        adjoint_error_estimate(plain, quadrature) rtol = 1.0e-6
    @test_throws ArgumentError adjoint_error_estimate(
        prob, GlobalAdjoint(Tsit5(); sensealg = InterpolatingAdjoint(checkpointing = true))
    )
    full = ODEProblem(
        ODEFunction(decay; mass_matrix = [1.0 0.5; 0.0 0.0]), [1.0, 1.0], (0.0, 1.0), [1.0]
    )
    @test_throws ArgumentError adjoint_error_estimate(full, GlobalAdjoint(Rodas5P()))
    # A scaled identity fails, because the error identity assumes M = I.
    scaled = ODEProblem(ODEFunction(decay; mass_matrix = 2.0I), [1.0, 1.0], (0.0, 1.0), [1.0])
    @test_throws ArgumentError adjoint_error_estimate(scaled, GlobalAdjoint(Rodas5P()))
    @test_throws ArgumentError solve(scaled, GlobalAdjoint(Rodas5P(); gtol = 1.0e-4))
end
