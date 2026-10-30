# The DefectCorrection companion with support for callbacks and mass matrices (for
# example ModelingToolkit models with events), in both modes. The companion is the true
# problem in the coordinates e = u - P(t):
#
#     M e' = f(P(t) + e, p̃, t) - M P'(t),    e(t₀) = 0,
#
# where P is the dense output of the forward solve and p̃ is the companion's own copy of
# the initial parameters. Then P + e is the exact solution of the true problem, so the
# companion needs three rules besides its equation:
#
#   - P and P' come from the curve of the current interval of P (see `_IntervalInterp`).
#   - P jumps at its knots, while the true solution does not jump there, so at each
#     knot e ← e - [P]. At an event knot the jump includes the forward affect; the true
#     solution gets its own affect from the user callbacks below.
#   - The user callbacks act on the true state P + e and on p̃: conditions are
#     evaluated there and affects are applied there, through `_TrueStateIntegrator`.
#     The true event times can differ from the forward ones; the companion follows the
#     true ones, and e is O(1) in the short gap between the two.
#
# A DiscreteCallback with a state condition is checked only at the knots of P, where the
# forward solve checked it, so that the true solution has a definition that does not
# depend on the step grid of the companion.
#
# InterpolatingMode integrates the companion after the forward solve, against its stored
# dense output (`_IntervalInterp`). SimultaneousMode advances one persistent companion
# integrator over each forward step, against the live interpolant of that step
# (`_LiveInterp`). Both give P(t), P'(t), whether t is a forward step end, and the
# forward solution, so the callback wrappers below serve both.

# P and P' on the current interval [ends[k - 1], ends[k]] of the forward solution: the
# right limit at the start of the interval and the left limit at its end. A
# right-continuous P is wrong at the end of a step that ends at an event knot: there it
# is the curve after the forward affect, while e still has its value from before the
# knot update, so the true state P + e would jump by the affect and events are missed.
# The knot callback advances `k`.
mutable struct _IntervalInterp{S, T}
    sol::S
    ends::Vector{T}
    k::Int
end

function (P::_IntervalInterp)(t)
    return t >= P.ends[P.k] ? P.sol(t; continuity = :left) : P.sol(t; continuity = :right)
end

function _interval_deriv(P::_IntervalInterp, t)
    return t >= P.ends[P.k] ? P.sol(t, Val{1}; continuity = :left) :
        P.sol(t, Val{1}; continuity = :right)
end

# In-place forms for the companion RHS, which runs at every stage.
_interval_continuity(P::_IntervalInterp, t) = t >= P.ends[P.k] ? :left : :right
function _interval_value!(out, P::_IntervalInterp, t)
    return P.sol(out, t; continuity = _interval_continuity(P, t))
end
function _interval_deriv!(out, P::_IntervalInterp, t)
    return P.sol(out, t, Val{1}; continuity = _interval_continuity(P, t))
end

_is_step_end(P::_IntervalInterp, t) = insorted(t, P.ends)
_forward_solution(P::_IntervalInterp) = P.sol

# P and P' on the current step of the forward integrator in SimultaneousMode.
struct _LiveInterp{I}
    integrator::I
end

(P::_LiveInterp)(t) = P.integrator(t)
_interval_deriv(P::_LiveInterp, t) = P.integrator(t, Val{1})
_interval_value!(out, P::_LiveInterp, t) = P.integrator(out, t)
_interval_deriv!(out, P::_LiveInterp, t) = P.integrator(out, t, Val{1})
_is_step_end(P::_LiveInterp, t) = t == P.integrator.t
_forward_solution(P::_LiveInterp) = P.integrator.sol

# The interior knots (unique saved times strictly inside the span) and the full jump of
# the dense output at each of them.
function _knot_jumps(sol)
    times = unique(sol.t)
    knots = times[2:(end - 1)]
    jumps = [sol(t; continuity = :right) .- sol(t; continuity = :left) for t in knots]
    return knots, jumps
end

function _knot_update!(integrator, P, knots, jumps)
    j = searchsortedfirst(knots, integrator.t)
    P.k == j || throw(
        ErrorException(
            "the companion solve lost track of the forward intervals at t = $(integrator.t)"
        )
    )
    integrator.u .-= jumps[j]
    P.k += 1
    return nothing
end

# The error trajectory. At each step end of the forward solve, the companion records the
# true state P(t) + e: before the user DiscreteCallbacks run there (after the knot update),
# after each of their calls, and after all of them. The true solution is continuous at a
# knot without a true event, so the values are equal there; a true event of a
# DiscreteCallback at the knot makes them different.
function _true_state_recorder(P, records, first::Bool)
    affect! = function (integrator)
        _record_true_state!(records, P, integrator, first)
        SciMLBase.u_modified!(integrator, false)
        return nothing
    end
    return SciMLBase.DiscreteCallback(
        (e, t, integrator) -> _is_step_end(P, t), affect!;
        save_positions = (false, false), initializealg = SciMLBase.NoInit()
    )
end

function _record_true_state!(records, P, integrator, first = false)
    t = integrator.t
    true_state = P(t) .+ integrator.u
    first ? get!(() -> [true_state], records, t) : push!(records[t], true_state)
    return nothing
end

# The recorders, or nothing without a trajectory: the first one runs after the knot update
# and before the user callbacks, the hook after each companion DiscreteCallback call, and
# the last recorder after all of them.
function _trajectory_recorders(P, trajectory, T)
    trajectory || return nothing, (), (), nothing
    records = Dict{T, Any}()
    return records, (_true_state_recorder(P, records, true),),
        (_true_state_recorder(P, records, false),),
        integrator -> _record_true_state!(records, P, integrator)
end

# The global error y(tᵢ) - uᵢ of each saved state of the forward solution. A callback can
# save before and after its affect, so an event time can be saved several times; `calls`
# (the forward calls with the number of saves before each one) places each entry between
# the calls at its time. An entry after k DiscreteCallback calls gets the true state after
# k companion calls there (the true counterpart of a continuous event is at another time,
# so the true state at t includes it or not). At t₀ the true solution starts at the saved
# state.
function _trajectory_errors(sol, records, calls)
    t = sol.t
    return map(eachindex(t)) do i
        t[i] == t[begin] && return zero(sol.u[i])
        discrete = [c for c in calls if c.t == t[i] && c.discrete]
        k = count(c -> c.saves < i, discrete)
        states = records[t[i]]
        true_state = k == length(discrete) && k > 0 ? states[end] : states[min(1 + k, end)]
        return true_state .- sol.u[i]
    end
end

# The forward calls, each with its time, whether it is a DiscreteCallback, and the number of
# saved states before it. `marks` holds (t, number of events before, number of saves before).
function _forward_call_saves(marks, events)
    calls = _forward_calls([(t, n) for (t, n, _) in marks], events)
    return [
        (; t, discrete = last(first(labels)) == 0, saves = marks[k][3])
            for (k, (t, labels)) in enumerate(calls)
    ]
end

# The view of the companion integrator that the user callbacks see: the true state
# u = P + e, the companion parameters p̃, and the time. Changes to `u` are written back
# to e by `_run_on_true_state!`. Integrator functions that have no meaning for the true
# state throw, so that an unsupported callback fails with a clear message.
mutable struct _TrueStateIntegrator{U, P, T, I, S, C}
    u::U
    p::P
    t::T
    parent::I
    sol::S
    interp::C
end

# The true state at a time in the current forward interval (a SavingCallback with `saveat`
# interpolates it): the curve of P on that interval (the left limit at its end) plus the
# companion interpolant of e. The interpolant covers only the last companion step; before
# it, it extrapolates, and the value is less accurate (about the size of the forward
# error). The companion runs on a copy of the user callbacks, so such values do not reach
# the user; only an affect that makes a decision from an interpolated state would see them.
_step_curve(P::_IntervalInterp, t) = P.sol(t; continuity = :left)
_step_curve(P::_LiveInterp, t) = P.integrator(t)
function (integrator::_TrueStateIntegrator)(t)
    return _step_curve(integrator.interp, t) .+ integrator.parent(t)
end
(integrator::_TrueStateIntegrator)(out, t) = (out .= integrator(t))

# Read-only properties of the time axis and the options mean the same for the true state,
# so they come from the companion integrator: ModelingToolkit's implicit affects read the
# tolerances, and PresetTimeCallback, PeriodicCallback, and IterativeCallback read `tdir`.
const _COMPANION_PROPERTIES = (:opts, :tdir, :dt, :tprev, :iter)

function Base.getproperty(integrator::_TrueStateIntegrator, name::Symbol)
    name in _COMPANION_PROPERTIES && return getproperty(getfield(integrator, :parent), name)
    return getfield(integrator, name)
end

function Base.propertynames(integrator::_TrueStateIntegrator, private::Bool = false)
    return (fieldnames(typeof(integrator))..., _COMPANION_PROPERTIES...)
end

SymbolicIndexingInterface.state_values(integrator::_TrueStateIntegrator) = integrator.u
SymbolicIndexingInterface.parameter_values(integrator::_TrueStateIntegrator) = integrator.p
SymbolicIndexingInterface.current_time(integrator::_TrueStateIntegrator) = integrator.t
SymbolicIndexingInterface.is_time_dependent(::_TrueStateIntegrator) = true
function SymbolicIndexingInterface.symbolic_container(integrator::_TrueStateIntegrator)
    return SymbolicIndexingInterface.symbolic_container(integrator.sol)
end
function SymbolicIndexingInterface.set_state!(integrator::_TrueStateIntegrator, val, idx)
    return integrator.u[idx] = val
end
function SymbolicIndexingInterface.set_parameter!(integrator::_TrueStateIntegrator, val, idx)
    return SymbolicIndexingInterface.set_parameter!(integrator.p, val, idx)
end

# Compiled ModelingToolkit conditions read `get_sol(integrator).prob`.
SciMLBase.get_sol(integrator::_TrueStateIntegrator) = integrator.sol
SciMLBase.get_tmp_cache(integrator::_TrueStateIntegrator) = (similar(integrator.u),)
SciMLBase.u_modified!(::_TrueStateIntegrator, ::Bool) = nothing
SciMLBase.derivative_discontinuity!(::_TrueStateIntegrator, ::Bool) = nothing
SciMLBase.reeval_internals_due_to_modification!(::_TrueStateIntegrator, args...; kwargs...) = nothing
# A tstop is a time, so it means the same for the true state (PresetTimeCallback).
function SciMLBase.add_tstop!(integrator::_TrueStateIntegrator, t)
    return SciMLBase.add_tstop!(integrator.parent, t)
end
# PeriodicCallback reads the last tstop before it adds the next one.
DiffEqBase.get_tstops_max(integrator::_TrueStateIntegrator) = DiffEqBase.get_tstops_max(integrator.parent)
# A terminating affect ends the forward solve at its event; the estimate is the error at that
# time, so the true state goes on to it (the companion ends with the forward solve).
SciMLBase.terminate!(::_TrueStateIntegrator, args...; kwargs...) = nothing
# The derivative of the true state of an ODE is f(u, p, t) (TerminateSteadyState reads it).
# For a DAE it has no value in the algebraic components.
function SciMLBase.get_du(integrator::_TrueStateIntegrator)
    prob = integrator.sol.prob
    prob.f.mass_matrix isa LinearAlgebra.UniformScaling || throw(
        ArgumentError(
            "GlobalErrorEstimation applies the user callbacks to the true state in its " *
                "companion solve; a callback that calls `get_du` is not supported for a " *
                "problem with a mass matrix"
        )
    )
    return _rhs(SciMLBase.unwrapped_f(prob.f), SciMLBase.isinplace(prob), integrator.u, integrator.p, integrator.t)
end
SciMLBase.get_du!(out, integrator::_TrueStateIntegrator) = (out .= SciMLBase.get_du(integrator))

for name in (
        :add_saveat!, :set_proposed_dt!, :get_proposed_dt, :reinit!,
        :auto_dt_reset!, :savevalues!, :set_u!, :set_t!,
    )
    @eval function SciMLBase.$name(::_TrueStateIntegrator, args...; kwargs...)
        throw(
            ArgumentError(
                "GlobalErrorEstimation applies the user callbacks to the true state in its " *
                    "companion solve; a callback that calls `" * $(string(name)) *
                    "` is not supported"
            )
        )
    end
end

function _true_state_integrator(ũ, integrator, P)
    return _TrueStateIntegrator(
        ũ, integrator.p, integrator.t, integrator, _forward_solution(P), P
    )
end

# Run `f!(true_state_integrator, args...)` on u = P(t) + e and write the result back to e.
function _run_on_true_state!(f!, integrator, P, args...)
    Pt = P(integrator.t)
    true_state = _true_state_integrator(Pt .+ integrator.u, integrator, P)
    f!(true_state, args...)
    integrator.u .= true_state.u .- Pt
    true_state.p === integrator.p || (integrator.p = true_state.p)
    return nothing
end

function _true_state_initialize(callback, P)
    return (_, e, t, integrator) -> _run_on_true_state!(
        true_state -> callback.initialize(callback, true_state.u, t, true_state),
        integrator, P
    )
end

# Each event is recorded as a label (callback number, kind), so that the forward and the
# companion event sequences can be compared. The kind is +1 or -1 for the up- and
# down-crossing affect of a ContinuousCallback, 0 for a DiscreteCallback, and the fired
# component for a VectorContinuousCallback, whose affect gets either the index of the
# fired component or a mask of all the components that fired at the same time.
_fired_components(index::Integer) = (index,)
_fired_components(mask) = findall(!iszero, mask)

function _record_event!(events, label, args...)
    callback_number, kind = label
    if kind === nothing
        for component in _fired_components(only(args))
            push!(events, (callback_number, component))
        end
    else
        push!(events, (callback_number, kind))
    end
    return nothing
end

# After an affect, the forward solve reinitializes a DAE with the initialization algorithm
# of the callback, or with the one of the solve if the callback has none. The companion
# applies the same reinitialization to the true state. The companion is a DAE in e with
# the algebraic equations 0 = g(P + e), so an algorithm that solves the equations of its
# integrator, such as BrownFullBasicInit, gives the reinitialized true state; its
# tolerance must be below the size of e. After CheckInit and NoInit the forward state is
# consistent without a change, and so is the true state. OverrideInit solves the
# ModelingToolkit initialization system for u, which has no meaning for e.
function _companion_reinit(alg, prob, kwargs, abstol)
    if alg === nothing
        problem_kwargs = values(prob.kwargs)
        alg = get(kwargs, :initializealg, get(problem_kwargs, :initializealg, DefaultInit()))
    end
    if alg isa DefaultInit
        alg = SciMLBase.has_initializeprob(prob.f) ? SciMLBase.OverrideInit() :
            SciMLBase.CheckInit()
    end
    return _companion_reinit(alg, abstol)
end
_companion_reinit(::Union{SciMLBase.NoInit, SciMLBase.CheckInit}, _) = SciMLBase.NoInit()
function _companion_reinit(alg::BrownFullBasicInit, abstol)
    return BrownFullBasicInit(min(alg.abstol, minimum(abstol)), alg.nlsolve)
end
_companion_reinit(alg, _) = alg

# The reinitialization runs in the affect wrapper, not as the initialization algorithm of
# the companion callback, because only there is it known whether the companion is a DAE.
_reinitialize_true_state!(integrator, ::SciMLBase.NoInit) = nothing
function _reinitialize_true_state!(integrator, alg)
    integrator.isdae || return nothing
    alg isa SciMLBase.OverrideInit && throw(
        ArgumentError(
            "the reinitialization of a DAE with OverrideInit after an event is not " *
                "supported in the companion error solve (OverrideInit is the default for " *
                "a ModelingToolkit problem when the callback has no `initializealg`). " *
                "Give the callback `initializealg = BrownFullBasicInit()`, `CheckInit()` " *
                "or `NoInit()`."
        )
    )
    SciMLBase.initialize_dae!(integrator, alg)
    integrator.sol.retcode == SciMLBase.ReturnCode.InitialFailure && throw(
        ErrorException("the companion reinitialization of the true state after an event failed")
    )
    return nothing
end

# The condition and the affects of a user callback, applied to the true state. The
# companion callbacks save nothing (the companion saves only its final value);
# `after_call(integrator)` runs after each call (the error trajectory records the true state).
function _recording_true_state_affect(
        affect!, P, events, label, reinit, after_call = nothing
    )
    affect! === nothing && return nothing
    return function (integrator, args...)
        _record_event!(events, label, args...)
        _run_on_true_state!(affect!, integrator, P, args...)
        _reinitialize_true_state!(integrator, reinit)
        after_call === nothing || after_call(integrator)
        return nothing
    end
end

# The companion callbacks have no finalizer: the user finalizers run in the forward solve.
_finalize_nothing(_, _, _, _) = nothing

_true_state_view(ũ, ::Nothing) = ũ
_true_state_view(ũ, idxs) = view(ũ, idxs)

function _companion_callback(
        callback::SciMLBase.ContinuousCallback, P, events, number, reinit, _ = nothing
    )
    condition = function (e, t, integrator)
        ũ = P(t) .+ e
        return callback.condition(
            _true_state_view(ũ, callback.idxs), t, _true_state_integrator(ũ, integrator, P)
        )
    end
    return SciMLBase.ContinuousCallback(
        condition,
        _recording_true_state_affect(callback.affect!, P, events, (number, 1), reinit),
        _recording_true_state_affect(callback.affect_neg!, P, events, (number, -1), reinit),
        _true_state_initialize(callback, P), _finalize_nothing, nothing,
        callback.rootfind, callback.interp_points, (false, false), callback.dtrelax,
        callback.abstol, callback.reltol, callback.repeat_nudge, SciMLBase.NoInit(),
        (), callback.maybe_discontinuity, false
    )
end

function _companion_callback(
        callback::SciMLBase.VectorContinuousCallback, P, events, number, reinit, _ = nothing
    )
    condition = function (out, e, t, integrator)
        ũ = P(t) .+ e
        return callback.condition(
            out, _true_state_view(ũ, callback.idxs), t,
            _true_state_integrator(ũ, integrator, P)
        )
    end
    return SciMLBase.VectorContinuousCallback(
        condition,
        _recording_true_state_affect(callback.affect!, P, events, (number, nothing), reinit),
        callback.len,
        _true_state_initialize(callback, P), _finalize_nothing, nothing,
        callback.rootfind, callback.interp_points, (false, false), callback.dtrelax,
        callback.abstol, callback.reltol, callback.repeat_nudge, SciMLBase.NoInit(),
        (), callback.maybe_discontinuity, false
    )
end

# The forward solve checks a DiscreteCallback only at its own step ends, so the
# companion checks it only at the knots of P and at the final time, not at the step ends
# of the companion. P(t) + e is the true state there on either side of the knot update.
function _companion_callback(
        callback::SciMLBase.DiscreteCallback, P, events, number, reinit, after_call = nothing
    )
    condition = function (e, t, integrator)
        _is_step_end(P, t) || return false
        ũ = P(t) .+ e
        return callback.condition(ũ, t, _true_state_integrator(ũ, integrator, P))
    end
    return SciMLBase.DiscreteCallback(
        condition,
        _recording_true_state_affect(callback.affect!, P, events, (number, 0), reinit, after_call),
        _true_state_initialize(callback, P), _finalize_nothing, (false, false),
        SciMLBase.NoInit(), (), false
    )
end

# The forward solve runs the user callbacks unchanged, except that each affect also
# records its event label, for the comparison with the companion events, and first calls
# `before(integrator)` (SimultaneousMode advances the companion there, while the
# interpolant of the step that the event cuts is still valid). `before` gets the label and
# the arguments of the affect too; if it returns a function, that runs after the affect.
function _recording_affect(affect!, events, label, before)
    affect! === nothing && return nothing
    return function (integrator, args...)
        after = before === nothing ? nothing : before(integrator, label, args...)
        _record_event!(events, label, args...)
        result = affect!(integrator, args...)
        after === nothing || after(integrator)
        return result
    end
end

# The user initialize and finalize functions get the user callback, not the recording one:
# some of them read its fields (the SavingCallback initialize reads `cb.affect!`).
_original_initialize(callback) =
    (_, u, t, integrator) -> callback.initialize(callback, u, t, integrator)
_original_finalize(callback) =
    (_, u, t, integrator) -> callback.finalize(callback, u, t, integrator)

function _recording_callback(callback::SciMLBase.ContinuousCallback, events, number, before)
    return SciMLBase.ContinuousCallback(
        callback.condition, _recording_affect(callback.affect!, events, (number, 1), before),
        _recording_affect(callback.affect_neg!, events, (number, -1), before),
        _original_initialize(callback),
        _original_finalize(callback), callback.idxs, callback.rootfind, callback.interp_points,
        callback.save_positions, callback.dtrelax, callback.abstol, callback.reltol,
        callback.repeat_nudge, callback.initializealg, callback.saved_clock_partitions,
        callback.maybe_discontinuity, callback.initialize_save_discretes
    )
end

function _recording_callback(callback::SciMLBase.VectorContinuousCallback, events, number, before)
    return SciMLBase.VectorContinuousCallback(
        callback.condition,
        _recording_affect(callback.affect!, events, (number, nothing), before),
        callback.len,
        _original_initialize(callback), _original_finalize(callback), callback.idxs,
        callback.rootfind, callback.interp_points, callback.save_positions, callback.dtrelax,
        callback.abstol, callback.reltol, callback.repeat_nudge, callback.initializealg,
        callback.saved_clock_partitions, callback.maybe_discontinuity,
        callback.initialize_save_discretes
    )
end

function _recording_callback(callback::SciMLBase.DiscreteCallback, events, number, before)
    return SciMLBase.DiscreteCallback(
        callback.condition, _recording_affect(callback.affect!, events, (number, 0), before),
        _original_initialize(callback), _original_finalize(callback), callback.save_positions,
        callback.initializealg, callback.saved_clock_partitions,
        callback.initialize_save_discretes
    )
end

_callback_list(callbacks) = (callbacks.continuous_callbacks..., callbacks.discrete_callbacks...)

# The user callbacks from the problem and from the solve keywords, as the solve would
# merge them.
function _user_callbacks(prob, kwargs)
    problem_kwargs = values(prob.kwargs)
    from_problem = haskey(problem_kwargs, :callback) ? problem_kwargs.callback : nothing
    from_solve = haskey(kwargs, :callback) ? kwargs[:callback] : nothing
    return SciMLBase.CallbackSet(from_problem, from_solve)
end

_without_callback(kwargs) = (; (k => v for (k, v) in pairs(kwargs) if k !== :callback)...)
_with_copied_callback(kwargs) =
    haskey(kwargs, :callback) ? (; kwargs..., callback = deepcopy(kwargs[:callback])) : kwargs

# The forward problem and callbacks. The forward affects change the parameters in place,
# so the forward solve gets its own copy (the user's problem keeps its parameters). The
# callbacks are passed explicitly, so they are removed from the problem keywords to avoid
# a second merge.
function _forward_setup(prob, callback_list, forward_events, before)
    forward_prob = @set prob.p = deepcopy(prob.p)
    forward_prob = @set forward_prob.kwargs = _without_callback(prob.kwargs)
    forward_callbacks = SciMLBase.CallbackSet(
        (
            _recording_callback(c, forward_events, n, before)
                for (n, c) in enumerate(callback_list)
        )...
    )
    return forward_prob, forward_callbacks
end

_subtract_mass_times!(de, M::LinearAlgebra.UniformScaling, v, _) = (de .-= M.λ .* v)
function _subtract_mass_times!(de, M, v, buffer)
    eltype(v) === eltype(buffer) || return de .-= M * v
    LinearAlgebra.mul!(buffer, M, v)
    return de .-= buffer
end

# The companion problem M e' = f(P + e, p̃, t) - M P', with p̃ a copy of the initial
# parameters. P, P' and P + e go into buffers. An implicit companion solver can pass dual
# numbers for e (the Jacobian) or for t (the time derivative of a Rosenbrock method);
# then new arrays of their type are used.
function _companion_problem(prob, P, tspan)
    f = SciMLBase.unwrapped_f(prob.f)
    M = prob.f.mass_matrix
    iip = SciMLBase.isinplace(prob)
    Pt, dP, ũ, MdP = (similar(prob.u0) for _ in 1:4)
    rhs! = function (de, e, p̃, t)
        if t isa eltype(Pt)
            _interval_value!(Pt, P, t)
            _interval_deriv!(dP, P, t)
            Pv, dPv = Pt, dP
        else
            Pv = P(t)
            dPv = _interval_deriv(P, t)
        end
        x = eltype(e) === eltype(ũ) === eltype(Pv) ? (ũ .= Pv .+ e) : Pv .+ e
        iip ? f(de, x, p̃, t) : (de .= f(x, p̃, t))
        _subtract_mass_times!(de, M, dPv, MdP)
        return nothing
    end
    return SciMLBase.ODEProblem{true, SciMLBase.FullSpecialize}(
        SciMLBase.ODEFunction{true, SciMLBase.FullSpecialize}(rhs!; mass_matrix = M),
        zero(prob.u0), tspan, deepcopy(prob.p)
    )
end

# The companion runs on its own copy of the user callbacks, made before the forward solve
# changes them, as it runs on its own copy of the parameters: it solves the problem a
# second time, so the state inside a callback (a counter, the values of a SavingCallback)
# starts again from the beginning, and the user's callback objects see only the forward
# solve. `reinit` gives the companion reinitialization for the initialization algorithm
# of a callback (see `_companion_reinit`); `after_call` goes to the DiscreteCallbacks.
function _companion_callbacks(
        callback_list, P, companion_events, reinit, extra...; after_call = nothing
    )
    return SciMLBase.CallbackSet(
        extra...,
        (
            _companion_callback(c, P, companion_events, n, reinit(c.initializealg), after_call)
                for (n, c) in enumerate(callback_list)
        )...
    )
end

function _check_companion_retcode(retcode)
    SciMLBase.successful_retcode(retcode) || throw(
        ErrorException("the companion error solve failed with retcode $retcode")
    )
    return nothing
end

function _warn_if_events_differ(name, forward_events, companion_events)
    sort(companion_events) == sort(forward_events) || @warn(
        "$name: the events of the forward solve and of the true solution, as the " *
            "companion computes it, are different (forward: $(length(forward_events)), " *
            "true: $(length(companion_events))). The estimate includes the difference, " *
            "but it is sensitive to the event times."
    )
    return nothing
end

# InterpolatingMode: dense forward solve, then the companion against its interpolant.
# Returns the 2-norm of the endpoint error estimate, or with `trajectory = true` the
# forward solution with the error of each saved state in `global_error`.
function _defect_correction_estimate(
        name, prob, inner_alg, companion_alg, args...;
        abstol, reltol, companion_abstol, companion_reltol, trajectory = false, kwargs...
    )
    _validate_estimation_problem(prob, name; events = true)
    callback_list = _callback_list(_user_callbacks(prob, kwargs))
    companion_list = deepcopy(callback_list)
    forward_events = Tuple{Int, Int}[]
    marks = Tuple{eltype(prob.tspan), Int, Int}[]
    record_call = trajectory ?
        (integrator, _...) -> _mark_call!(marks, forward_events, integrator) : nothing
    forward_prob, forward_callbacks = _forward_setup(prob, callback_list, forward_events, record_call)
    solve_kwargs = merge(_without_callback(kwargs), _DENSE_SOLVE_KWARGS)
    forward = SciMLBase.init(
        forward_prob, inner_alg, args...;
        abstol, reltol, callback = forward_callbacks, solve_kwargs...
    )
    u_init = copy(forward.u)
    sol = SciMLBase.solve!(forward)
    _validate_estimation_solution(sol, name)

    knots, jumps = _knot_jumps(sol)
    T = sol.t[end]
    P = _IntervalInterp(sol, [knots; T], 1)
    companion_events = Tuple{Int, Int}[]
    knot_callback = SciMLBase.DiscreteCallback(
        (e, t, integrator) -> insorted(t, knots),
        integrator -> _knot_update!(integrator, P, knots, jumps);
        save_positions = (false, false), initializealg = SciMLBase.NoInit()
    )
    records, first_recorder, last_recorder, after_call =
        _trajectory_recorders(P, trajectory, eltype(sol.t))
    companion = SciMLBase.init(
        _companion_problem(prob, P, (prob.tspan[1], T)), companion_alg;
        abstol = companion_abstol, reltol = companion_reltol,
        callback = SciMLBase.CallbackSet(
            _companion_callbacks(
                companion_list, P, companion_events,
                alg -> _companion_reinit(alg, prob, kwargs, companion_abstol),
                knot_callback, first_recorder...; after_call
            ),
            last_recorder...
        ),
        tstops = knots, dense = false, save_everystep = false, save_start = false,
        save_end = true, initializealg = SciMLBase.NoInit()
    )
    # The initializations of the callbacks are exact (the true solution starts at the same
    # state), so the error after them is the state of the forward integrator after its
    # initialization minus P(t₀). The companion's own initializations set the internal state
    # of its copies of the callbacks; a change that they make to e is replaced here.
    SciMLBase.set_u!(companion, u_init .- P(prob.tspan[1]))
    companion_sol = SciMLBase.solve!(companion)
    _check_companion_retcode(companion_sol.retcode)
    _warn_if_events_differ(name, forward_events, companion_events)

    # The companion gives the error against P(T⁻); the saved state u(T) can differ from
    # P(T⁻), for example after an event at the final time.
    endpoint_error = companion_sol.u[end] .- (sol.u[end] .- sol(T; continuity = :left))
    trajectory || return LinearAlgebra.norm(endpoint_error)
    calls = _forward_call_saves(marks, forward_events)
    return @set sol.global_error = _trajectory_errors(sol, records, calls)
end

# Record a forward call: its time, the number of events before it, and the number of saved
# states before it (see `_trajectory_errors`).
function _mark_call!(marks, events, integrator)
    push!(marks, (integrator.t, length(events), length(integrator.sol.t)))
    return nothing
end

# SimultaneousMode: the state of the streaming companion between forward steps.
mutable struct _StreamingCompanion{C, V}
    companion::C
    # The end value of the curve of the last completed forward step, P(t⁻), once the
    # first step is done.
    last_end::V
    started::Bool
    storage::IdDict{Any, Any}
end

# A callback can keep its state in task-local storage under a key that its copies share
# (the DiffEqCallbacks PeriodicCallback does), so in SimultaneousMode, where the forward solve
# and the companion run in one task, the companion's copies of the callbacks get their own
# entries. The companion initialization runs with an empty storage (`_with_task_storage`),
# which shows the keys that its callbacks use; each companion step then swaps only those
# keys (`_with_companion_keys`, nothing for most callbacks), so other task-local state stays
# visible.
function _with_task_storage(f, storage)
    tls = task_local_storage()
    outer = copy(tls)
    empty!(tls)
    merge!(tls, storage)
    try
        return f()
    finally
        empty!(storage)
        merge!(storage, tls)
        empty!(tls)
        merge!(tls, outer)
    end
end

const _NO_ENTRY = gensym(:no_entry)

function _with_companion_keys(f, storage)
    isempty(storage) && return f()
    tls = task_local_storage()
    outer = Dict{Any, Any}(k => get(tls, k, _NO_ENTRY) for k in keys(storage))
    merge!(tls, storage)
    try
        return f()
    finally
        for k in keys(storage)
            storage[k] = tls[k]
        end
        for (k, v) in outer
            v === _NO_ENTRY ? delete!(tls, k) : (tls[k] = v)
        end
    end
end

# Advance the companion over the forward step that just completed, [tprev, t], against
# the live interpolant of that step. First subtract the jump of P at tprev: the start of
# the new curve minus the end of the old one (at an event this includes the forward
# affect). The companion stops at t, so its DiscreteCallbacks are checked there.
function _advance_streaming_companion!(state, forward)
    companion = state.companion
    t = forward.t
    companion.t < t || return nothing
    companion.t == forward.tprev || throw(
        ErrorException(
            "the companion solve lost track of the forward steps at t = $(forward.tprev)"
        )
    )
    if state.started
        SciMLBase.set_u!(companion, companion.u .- (forward(forward.tprev) .- state.last_end))
    end
    _with_companion_keys(() -> SciMLBase.step!(companion, t - companion.t, true), state.storage)
    _check_companion_retcode(SciMLBase.check_error(companion))
    state.last_end = forward(t)
    state.started = true
    return nothing
end

function _defect_correction_estimate_streaming(
        name, prob, inner_alg, companion_alg, args...;
        abstol, reltol, companion_abstol, companion_reltol, trajectory = false, kwargs...
    )
    _validate_estimation_problem(prob, name; events = true)
    callback_list = _callback_list(_user_callbacks(prob, kwargs))
    companion_list = deepcopy(callback_list)
    forward_events = Tuple{Int, Int}[]
    companion_events = Tuple{Int, Int}[]
    state_ref = Ref{Any}(nothing)
    # A forward affect advances the companion to the event time first. An affect during
    # the forward initialization (before the companion exists) needs no advance.
    marks = Tuple{eltype(prob.tspan), Int, Int}[]
    before = function (forward, _...)
        trajectory && _mark_call!(marks, forward_events, forward)
        state_ref[] === nothing || _advance_streaming_companion!(state_ref[], forward)
        return nothing
    end
    forward_prob, forward_callbacks = _forward_setup(prob, callback_list, forward_events, before)
    # The error trajectory needs the saved state at each forward step.
    save_kwargs = trajectory ? merge(_without_callback(kwargs), _DENSE_SOLVE_KWARGS) :
        (;
            dense = true, save_everystep = false, save_start = true, save_end = true,
            _without_callback(kwargs)...,
        )
    forward = SciMLBase.init(
        forward_prob, inner_alg, args...;
        abstol, reltol, callback = forward_callbacks, save_kwargs...
    )
    _validate_streaming_interpolation(forward.sol, name)

    P = _LiveInterp(forward)
    records, first_recorder, last_recorder, after_call = _trajectory_recorders(
        P, trajectory, eltype(forward.sol.t)
    )
    storage = IdDict{Any, Any}()
    companion = _with_task_storage(storage) do
        SciMLBase.init(
            _companion_problem(prob, P, prob.tspan), companion_alg;
            abstol = companion_abstol, reltol = companion_reltol,
            callback = SciMLBase.CallbackSet(
                _companion_callbacks(
                    companion_list, P, companion_events,
                    alg -> _companion_reinit(alg, prob, kwargs, companion_abstol),
                    first_recorder...; after_call
                ),
                last_recorder...
            ), dense = false,
            save_everystep = false, save_start = false, save_end = false,
            initializealg = SciMLBase.NoInit()
        )
    end
    # As in InterpolatingMode: the error after the initializations is the forward state after
    # them minus P(t₀), here zero (the first step of the live P starts at that state).
    SciMLBase.set_u!(companion, zero(companion.u))
    state = _StreamingCompanion(companion, copy(prob.u0), false, storage)
    state_ref[] = state
    for _ in forward
        _advance_streaming_companion!(state, forward)
    end
    SciMLBase.successful_retcode(forward.sol) ||
        throw(ErrorException("the forward solve failed with retcode $(forward.sol.retcode)"))
    _warn_if_events_differ(name, forward_events, companion_events)

    # As in InterpolatingMode: the error against P(T⁻), corrected to the saved state.
    endpoint_error = companion.u .- (forward.u .- state.last_end)
    trajectory || return LinearAlgebra.norm(endpoint_error)
    sol = forward.sol
    calls = _forward_call_saves(marks, forward_events)
    return @set sol.global_error = _trajectory_errors(sol, records, calls)
end
