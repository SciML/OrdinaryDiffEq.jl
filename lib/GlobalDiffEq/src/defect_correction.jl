# The DefectCorrection companion with callbacks and mass matrices (for example
# ModelingToolkit models with events), in both modes. The companion is the true problem in
# the coordinates e = u - P(t):
#
#     M e' = f(P(t) + e, p̃, t) - M P'(t),    e(t₀) = 0,
#
# where P is the dense output of the forward solve and p̃ is a companion copy of the initial
# parameters. Then P + e is the exact solution of the true problem. Besides its equation,
# the companion needs three rules:
#
#   - P and P' come from the curve of the current interval of P (see `_IntervalInterp`).
#   - P jumps at its knots and the true solution does not, so e ← e - [P] at each knot.
#     At an event knot, [P] includes the forward affect; the user callbacks below give the
#     true solution its own affect.
#   - The user callbacks see the true state P + e and p̃ (`_TrueStateIntegrator`). The
#     companion follows the true event times, which can differ from the forward ones; e is
#     O(1) in the short gap between the two.
#
# A DiscreteCallback with a state condition is checked only at the knots of P, as in the
# forward solve, so the true solution does not depend on the companion step grid. It is
# also checked at the true event times, because a continuous event ends a step in every
# solve.
#
# InterpolatingMode integrates the companion after the forward solve, against the stored
# dense output (`_IntervalInterp`). SimultaneousMode advances one companion integrator over
# each forward step, against the live interpolant of that step (`_LiveInterp`). Both types
# give P(t), P'(t), the forward step ends and the forward solution, so the callback
# wrappers below serve both.

# P and P' on the current interval [ends[k - 1], ends[k]] of the forward solution: the
# right limit at the interval start and the left limit at its end. The left limit matters
# at a step end on an event knot. There e still has its value from before the knot update,
# so P after the forward affect would shift the true state by the affect, and the companion
# would miss events. The knot callback advances `k`.
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

# The interior knots (unique saved times strictly inside the span) and the jump [P] of the
# dense output at each one.
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

# The error trajectory. At each forward step end, the companion records P(t) + e after the
# knot update (before the user DiscreteCallbacks), after each of their calls, and after all
# of them. The values differ only if a DiscreteCallback has a true event at the knot.
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

# The trajectory recorders, or nothing without a trajectory. The first recorder runs before
# the user callbacks and the last one after them. `after_call` runs after each companion
# DiscreteCallback call; it records only at a knot, because only knots have entries.
function _trajectory_recorders(P, trajectory, T)
    trajectory || return nothing, (), (), nothing
    records = Dict{T, Any}()
    after_call = function (integrator)
        _is_step_end(P, integrator.t) && _record_true_state!(records, P, integrator)
        return nothing
    end
    return records, (_true_state_recorder(P, records, true),),
        (_true_state_recorder(P, records, false),), after_call
end

# The global error y(tᵢ) - uᵢ of each saved forward state. A callback can save before and
# after its affect, so an event time can have several saves; `calls` places each save
# between the calls at its time. A save after k DiscreteCallback calls gets the true state
# after k companion calls there. A continuous event needs no such match, because its true
# counterpart is at another time. At t₀ the error is zero.
function _trajectory_errors(sol, records, calls)
    t = sol.t
    # The number of saves before each DiscreteCallback call, grouped by time.
    discrete_saves = Dict{eltype(t), Vector{Int}}()
    for c in calls
        c.discrete && push!(get!(discrete_saves, c.t, Int[]), c.saves)
    end
    return map(eachindex(t)) do i
        t[i] == t[begin] && return zero(sol.u[i])
        discrete = get(discrete_saves, t[i], Int[])
        k = count(<(i), discrete)
        states = records[t[i]]
        true_state = k == length(discrete) && k > 0 ? states[end] : states[min(1 + k, end)]
        return true_state .- sol.u[i]
    end
end

# The forward calls, each with its time, whether it is a DiscreteCallback, and the number
# of saves before it. `marks` holds (t, events before, saves before).
function _forward_call_saves(marks, events)
    calls = _forward_calls([(t, n) for (t, n, _) in marks], events)
    return [
        (; t, discrete = _is_discrete(first(labels)), saves = marks[k][3])
            for (k, (t, labels)) in enumerate(calls)
    ]
end

# The integrator that the user callbacks see: the true state u = P + e, the companion
# parameters p̃, and the time. `_run_on_true_state!` writes changes to `u` back to e.
# Integrator functions without a meaning for the true state throw a clear error.
mutable struct _TrueStateIntegrator{U, P, T, I, S, C}
    u::U
    p::P
    t::T
    parent::I
    sol::S
    interp::C
    # Set by `terminate!` (see below).
    terminated::Bool
end
function _TrueStateIntegrator(u, p, t, parent, sol, interp)
    return _TrueStateIntegrator(u, p, t, parent, sol, interp, false)
end

# The true state at a time in the current forward interval (a SavingCallback with `saveat`
# asks for it): the curve of P on that interval plus the companion interpolant of e. That
# interpolant covers only the last companion step and extrapolates before it, so values
# there are off by about the forward error. The user never sees such values, because the
# companion runs on a copy of the user callbacks.
_step_curve(P::_IntervalInterp, t) = P.sol(t; continuity = :left)
_step_curve(P::_LiveInterp, t) = P.integrator(t)
function (integrator::_TrueStateIntegrator)(t)
    return _step_curve(integrator.interp, t) .+ integrator.parent(t)
end
(integrator::_TrueStateIntegrator)(out, t) = (out .= integrator(t))

# Read-only properties of the time axis and the options mean the same for the true state,
# so they come from the companion integrator. ModelingToolkit implicit affects read the
# tolerances; PresetTimeCallback, PeriodicCallback and IterativeCallback read `tdir`.
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
# A terminating affect ends the forward solve at its event, and the estimate is the error
# at that time. So the companion does not stop; it only records the call for
# `_warn_if_events_differ`.
function SciMLBase.terminate!(integrator::_TrueStateIntegrator, args...; kwargs...)
    integrator.terminated = true
    return nothing
end
# The derivative of the true state of an ODE is f(u, p, t) (TerminateSteadyState reads it).
# A DAE has no derivative in the algebraic components.
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
# Returns whether `f!` called `terminate!`.
function _run_on_true_state!(f!, integrator, P, args...)
    Pt = P(integrator.t)
    true_state = _true_state_integrator(Pt .+ integrator.u, integrator, P)
    f!(true_state, args...)
    integrator.u .= true_state.u .- Pt
    true_state.p === integrator.p || (integrator.p = true_state.p)
    return true_state.terminated
end

function _true_state_initialize(callback, P)
    return (_, e, t, integrator) -> _run_on_true_state!(
        true_state -> callback.initialize(callback, true_state.u, t, true_state),
        integrator, P
    )
end

# An event label is (callback number, kind), for the comparison of the forward and the
# companion events. The kind is +1 or -1 for the up- and down-crossing affect of a
# ContinuousCallback and 0 for a DiscreteCallback. For a VectorContinuousCallback it is the
# fired component; the affect gets one index or a mask of the components that fired
# together.
_fired_components(index::Integer) = (index,)
_fired_components(mask) = findall(!iszero, mask)

# A DiscreteCallback label has kind 0. A VectorContinuousCallback label has kind `nothing`
# until `_record_event!` expands it.
_is_discrete(label) = last(label) === 0

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

# The companion form of the DAE reinitialization that the forward solve applies after an
# affect (the callback algorithm, else the solve algorithm). The companion is a DAE in e
# with the algebraic equations 0 = g(P + e), so BrownFullBasicInit gives the reinitialized
# true state; its tolerance must be below the size of e. CheckInit and NoInit change
# neither the forward nor the true state. OverrideInit solves for u, which has no meaning
# for e.
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

# The reinitialization runs in the affect wrapper, not as the callback `initializealg`,
# because only the wrapper knows if the companion is a DAE.
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

# A user affect applied to the true state. `record` collects the companion events, the
# times of the terminating calls, and the time of the last continuous event. `after_call`
# runs after each call (the trajectory recorder).
function _recording_true_state_affect(
        affect!, P, record, label, reinit, after_call = nothing
    )
    affect! === nothing && return nothing
    return function (integrator, args...)
        _record_event!(record.events, label, args...)
        _run_on_true_state!(affect!, integrator, P, args...) &&
            push!(record.terminations, integrator.t)
        _is_discrete(label) || (record.event_time[] = integrator.t)
        _reinitialize_true_state!(integrator, reinit)
        after_call === nothing || after_call(integrator)
        return nothing
    end
end

# The companion callbacks have no finalizer: the user finalizers run in the forward solve.
_finalize_nothing(_, _, _, _) = nothing

# The fields that differ from the user callback. The initialize runs on the true state, the
# user finalizers run in the forward solve, and the companion saves only its final value.
# `reinit` reinitializes the true state, so the companion state needs no `initializealg`.
# The condition applies `idxs` itself. Root finding and tolerances stay the same.
function _companion_fields(callback, P)
    fields = (;
        initialize = _true_state_initialize(callback, P), finalize = _finalize_nothing,
        save_positions = (false, false), initializealg = SciMLBase.NoInit(),
        saved_clock_partitions = (), initialize_save_discretes = false,
    )
    # A DiscreteCallback has no `idxs` field.
    return callback isa SciMLBase.DiscreteCallback ? fields : (; idxs = nothing, fields...)
end

_true_state_view(ũ, ::Nothing) = ũ
_true_state_view(ũ, idxs) = view(ũ, idxs)

function _companion_callback(
        callback::SciMLBase.ContinuousCallback, P, record, number, reinit, _ = nothing
    )
    condition = function (e, t, integrator)
        ũ = P(t) .+ e
        return callback.condition(
            _true_state_view(ũ, callback.idxs), t, _true_state_integrator(ũ, integrator, P)
        )
    end
    return setproperties(
        callback, (;
            condition,
            affect! = _recording_true_state_affect(callback.affect!, P, record, (number, 1), reinit),
            affect_neg! = _recording_true_state_affect(callback.affect_neg!, P, record, (number, -1), reinit),
            _companion_fields(callback, P)...,
        )
    )
end

function _companion_callback(
        callback::SciMLBase.VectorContinuousCallback, P, record, number, reinit, _ = nothing
    )
    condition = function (out, e, t, integrator)
        ũ = P(t) .+ e
        return callback.condition(
            out, _true_state_view(ũ, callback.idxs), t,
            _true_state_integrator(ũ, integrator, P)
        )
    end
    return setproperties(
        callback, (;
            condition,
            affect! = _recording_true_state_affect(callback.affect!, P, record, (number, nothing), reinit),
            _companion_fields(callback, P)...,
        )
    )
end

# The forward solve checks a DiscreteCallback only at its own step ends, so the companion
# checks it only at the knots of P and at the final time. P(t) + e is the true state there
# on either side of the knot update. A continuous event ends a step in every solve, so the
# companion checks a DiscreteCallback at the true event time (`record.event_time`), not at
# the forward event time (`record.event_knots`).
function _companion_callback(
        callback::SciMLBase.DiscreteCallback, P, record, number, reinit, after_call = nothing
    )
    condition = function (e, t, integrator)
        knot = _is_step_end(P, t) && !(t in record.event_knots)
        knot || record.event_time[] == t || return false
        ũ = P(t) .+ e
        return callback.condition(ũ, t, _true_state_integrator(ũ, integrator, P))
    end
    return setproperties(
        callback, (;
            condition,
            affect! = _recording_true_state_affect(callback.affect!, P, record, (number, 0), reinit, after_call),
            _companion_fields(callback, P)...,
        )
    )
end

# The forward solve runs the user callbacks unchanged, except that each affect records its
# event label and first calls `before(integrator, label, args...)`. SimultaneousMode
# advances the companion there, while the interpolant of the cut step is still valid. If
# `before` returns a function, it runs after the affect.
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

# The initialize and finalize of a recording callback, which call the user functions.
_original_hooks(callback) =
    (; initialize = _original_initialize(callback), finalize = _original_finalize(callback))

function _recording_callback(callback::SciMLBase.ContinuousCallback, events, number, before)
    return setproperties(
        callback, (;
            affect! = _recording_affect(callback.affect!, events, (number, 1), before),
            affect_neg! = _recording_affect(callback.affect_neg!, events, (number, -1), before),
            _original_hooks(callback)...,
        )
    )
end

function _recording_callback(callback::SciMLBase.VectorContinuousCallback, events, number, before)
    return setproperties(
        callback, (;
            affect! = _recording_affect(callback.affect!, events, (number, nothing), before),
            _original_hooks(callback)...,
        )
    )
end

function _recording_callback(callback::SciMLBase.DiscreteCallback, events, number, before)
    return setproperties(
        callback, (;
            affect! = _recording_affect(callback.affect!, events, (number, 0), before),
            _original_hooks(callback)...,
        )
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

_with_copied_callback(kwargs) =
    haskey(kwargs, :callback) ? (; kwargs..., callback = deepcopy(kwargs[:callback])) : kwargs

# The forward problem and callbacks. The forward affects change the parameters in place, so
# the forward solve gets its own copy. The callbacks are passed explicitly, so they are
# removed from the problem keywords to avoid a second merge.
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
# parameters. P, P' and P + e go into buffers. An implicit solver can pass dual numbers for
# e (the Jacobian) or for t (the Rosenbrock time derivative); then new arrays of that type
# hold them.
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

# The companion solves the problem a second time, so it gets its own copy of the user
# callbacks, made before the forward solve changes them. The callback state (a counter, the
# values of a SavingCallback) then starts from the beginning, and the user's callbacks see
# only the forward solve. `reinit` maps a callback `initializealg` to the companion
# reinitialization (`_companion_reinit`). Only the DiscreteCallbacks use `after_call`.
function _companion_callbacks(
        callback_list, P, record, reinit, extra...; after_call = nothing
    )
    return SciMLBase.CallbackSet(
        extra...,
        (
            _companion_callback(c, P, record, n, reinit(c.initializealg), after_call)
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

# Warn if the forward events differ from the true events that the companion computes.
# `marks` holds the forward calls as (t, events before, saves before). If the forward solve
# terminates, the companion can miss its continuous events at the end time, because the
# true event can come later. The same holds for the DiscreteCallback events that they
# trigger. Other DiscreteCallback events at the end time must agree. Also warn if only one
# of the two solves terminates at an event.
function _warn_if_events_differ(name, forward_events, record, marks, retcode, t_end, verbose)
    companion_events = record.events
    terminated = retcode == SciMLBase.ReturnCode.Terminated
    final_call = findfirst(m -> m[1] == t_end, marks)
    first_final = final_call === nothing ? length(forward_events) + 1 : marks[final_call][2] + 1
    final = forward_events[first_final:end]
    first_continuous = something(findfirst(!_is_discrete, final), length(final) + 1)
    required = terminated ?
        [forward_events[1:(first_final - 1)]; final[1:(first_continuous - 1)]] : forward_events
    same = sort(companion_events) == sort(forward_events)
    same || sort(companion_events) == sort(required) || @SciMLMessage(
        "$name: the events of the forward solve and of the true solution, as the " *
            "companion computes it, are different (forward: $(length(forward_events)), " *
            "true: $(length(companion_events))). The estimate includes the difference, " *
            "but it is sensitive to the event times.",
        verbose, :global_error_event_mismatch
    )
    if !terminated && !isempty(record.terminations)
        @SciMLMessage(
            "$name: the true solution, as the companion computes it, terminates at " *
                "t = $(first(record.terminations)), but the forward solve does not. The " *
                "estimate is the error at the end time of the forward solve.",
            verbose, :global_error_event_mismatch
        )
    elseif terminated && same && isempty(record.terminations)
        @SciMLMessage(
            "$name: the forward solve terminates at t = $t_end, but the true solution, as " *
                "the companion computes it, does not terminate at the same event. The " *
                "estimate is the error at the end time of the forward solve.",
            verbose, :global_error_event_mismatch
        )
    end
    return nothing
end

# The forward part of a DefectCorrection estimate, shared by both modes. The companion copy
# of the callbacks is made before the forward solve changes them. `hook(marks, events)`
# gives the mode-specific rest of the `before` function (`_forward_setup`). The times are
# floats, also for an integer time span.
function _defect_correction_forward(name, prob, kwargs, hook)
    _validate_estimation_problem(prob, name; events = true)
    callback_list = _callback_list(_user_callbacks(prob, kwargs))
    companion_list = deepcopy(callback_list)
    events = Tuple{Int, Int}[]
    T = float(eltype(prob.tspan))
    marks = Tuple{T, Int, Int}[]
    event_knots = Set{T}()
    rest = hook(marks, events)
    before = function (integrator, label, args...)
        _is_discrete(label) || push!(event_knots, integrator.t)
        return rest(integrator, label, args...)
    end
    forward_prob, forward_callbacks = _forward_setup(prob, callback_list, events, before)
    return (; companion_list, events, marks, event_knots, forward_prob, forward_callbacks)
end

# The companion callbacks: the user callbacks on the true state, the `extra` callbacks (the
# InterpolatingMode knot update) and, with `trajectory`, the true-state recorders. `record`
# also holds the forward continuous event times (see `_companion_callback`).
function _defect_correction_callbacks(forward, P, prob, kwargs, companion_abstol, trajectory, extra...)
    T = float(eltype(prob.tspan))
    record = (;
        events = Tuple{Int, Int}[], terminations = T[], event_time = Ref(T(NaN)),
        forward.event_knots,
    )
    records, first_recorder, last_recorder, after_call = _trajectory_recorders(P, trajectory, T)
    callback = SciMLBase.CallbackSet(
        _companion_callbacks(
            forward.companion_list, P, record,
            alg -> _companion_reinit(alg, prob, kwargs, companion_abstol),
            extra..., first_recorder...; after_call
        ),
        last_recorder...
    )
    return (; record, records, callback)
end

# The 2-norm of the endpoint error, or with `trajectory` also the forward solution with the
# error of each saved state in `global_error`, as `(; estimate, sol)`. It also warns about
# event differences. `t_end` is the forward end time; a solution that saves only some steps
# can end before it.
function _defect_correction_result(
        name, forward, callbacks, sol, endpoint_error, trajectory; verbose, t_end = sol.t[end]
    )
    _warn_if_events_differ(
        name, forward.events, callbacks.record, forward.marks, sol.retcode, t_end, verbose
    )
    estimate = LinearAlgebra.norm(endpoint_error)
    trajectory || return estimate
    calls = _forward_call_saves(forward.marks, forward.events)
    return (; estimate, sol = @set sol.global_error = _trajectory_errors(sol, callbacks.records, calls))
end

# InterpolatingMode: a dense forward solve, then the companion against its interpolant.
# The result is that of `_defect_correction_result`.
function _defect_correction_estimate(
        name, prob, inner_alg, companion_alg, args...;
        abstol, reltol, companion_abstol, companion_reltol, trajectory = false, kwargs...
    )
    forward = _defect_correction_forward(
        name, prob, kwargs,
        (marks, events) -> (integrator, _...) -> _mark_call!(marks, events, integrator)
    )
    integrator = SciMLBase.init(
        forward.forward_prob, inner_alg, args...;
        abstol, reltol, callback = forward.forward_callbacks,
        merge(_without_callback(kwargs), _DENSE_SOLVE_KWARGS)...
    )
    u_init = copy(integrator.u)
    sol = SciMLBase.solve!(integrator)
    _validate_estimation_solution(sol, name)

    knots, jumps = _knot_jumps(sol)
    T = sol.t[end]
    P = _IntervalInterp(sol, [knots; T], 1)
    knot_callback = SciMLBase.DiscreteCallback(
        (e, t, integrator) -> insorted(t, knots),
        integrator -> _knot_update!(integrator, P, knots, jumps);
        save_positions = (false, false), initializealg = SciMLBase.NoInit()
    )
    callbacks = _defect_correction_callbacks(
        forward, P, prob, kwargs, companion_abstol, trajectory, knot_callback
    )
    companion = SciMLBase.init(
        _companion_problem(prob, P, (prob.tspan[1], T)), companion_alg;
        abstol = companion_abstol, reltol = companion_reltol, callbacks.callback,
        tstops = knots, dense = false, save_everystep = false, save_start = false,
        save_end = true, initializealg = SciMLBase.NoInit()
    )
    # The callback initializations are exact (the true solution starts at the same state),
    # so e after them is the initialized forward state minus P(t₀). The companion
    # initializations only set the state inside the callback copies; `set_u!` replaces any
    # change that they make to e.
    SciMLBase.set_u!(companion, u_init .- P(prob.tspan[1]))
    companion_sol = SciMLBase.solve!(companion)
    _check_companion_retcode(companion_sol.retcode)

    # The companion gives the error against P(T⁻); the saved state u(T) can differ from
    # P(T⁻), for example after an event at the final time.
    endpoint_error = companion_sol.u[end] .- (sol.u[end] .- sol(T; continuity = :left))
    return _defect_correction_result(
        name, forward, callbacks, sol, endpoint_error, trajectory; verbose = _verbosity(kwargs)
    )
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
    # P(t⁻) at the end of the last completed forward step, after the first step.
    last_end::V
    started::Bool
    storage::IdDict{Any, Any}
end

# Some callbacks keep state in task-local storage under a key that their copies share (the
# DiffEqCallbacks PeriodicCallback does). In SimultaneousMode both solves run in one task,
# so the companion copies need their own entries. The companion initialization runs with an
# empty storage, which shows the keys of its callbacks (`_with_task_storage`). Each
# companion step swaps only those keys (`_with_companion_keys`), so other task-local state
# stays visible.
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

# Advance the companion over the forward step [tprev, t] that just completed, against the
# live interpolant of that step. First subtract the jump of P at tprev; at an event it
# includes the forward affect. The companion stops at t, so its DiscreteCallbacks are
# checked there.
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
    _with_companion_keys(() -> _step_to!(companion, t), state.storage)
    _check_companion_retcode(SciMLBase.check_error(companion))
    state.last_end = forward(t)
    state.started = true
    return nothing
end

# Step the integrator to exactly t, the forward step end, with a tstop at t.
# `step!(integrator, t - integrator.t, true)` can miss t by rounding if t > 2 integrator.t.
function _step_to!(integrator, t)
    SciMLBase.add_tstop!(integrator, t)
    while integrator.t < t
        SciMLBase.step!(integrator)
        integrator.sol.retcode in (SciMLBase.ReturnCode.Default, SciMLBase.ReturnCode.Success) ||
            break
    end
    return nothing
end

function _defect_correction_estimate_streaming(
        name, prob, inner_alg, companion_alg, args...;
        abstol, reltol, companion_abstol, companion_reltol, trajectory = false, kwargs...
    )
    # A forward affect first advances the companion to the event time. An affect in the
    # forward initialization needs no advance, because the companion does not exist yet.
    state_ref = Ref{Any}(nothing)
    forward = _defect_correction_forward(
        name, prob, kwargs,
        (marks, events) -> function (integrator, _...)
            # Without a trajectory, `_warn_if_events_differ` needs only the marks at the
            # latest call time, so the memory stays O(1).
            trajectory || isempty(marks) || first(last(marks)) == integrator.t || empty!(marks)
            _mark_call!(marks, events, integrator)
            state_ref[] === nothing || _advance_streaming_companion!(state_ref[], integrator)
            return nothing
        end
    )
    # The error trajectory needs the saved state at each forward step.
    save_kwargs = trajectory ? merge(_without_callback(kwargs), _DENSE_SOLVE_KWARGS) :
        (;
            dense = true, save_everystep = false, save_start = true, save_end = true,
            _without_callback(kwargs)...,
        )
    integrator = SciMLBase.init(
        forward.forward_prob, inner_alg, args...;
        abstol, reltol, callback = forward.forward_callbacks, save_kwargs...
    )
    _validate_streaming_interpolation(integrator.sol, name)

    P = _LiveInterp(integrator)
    callbacks = _defect_correction_callbacks(forward, P, prob, kwargs, companion_abstol, trajectory)
    storage = IdDict{Any, Any}()
    companion = _with_task_storage(storage) do
        SciMLBase.init(
            _companion_problem(prob, P, prob.tspan), companion_alg;
            abstol = companion_abstol, reltol = companion_reltol, callbacks.callback,
            dense = false, save_everystep = false, save_start = false, save_end = false,
            initializealg = SciMLBase.NoInit()
        )
    end
    # As in InterpolatingMode, e after the initializations is the forward state minus
    # P(t₀). That is zero here, because the first step of the live P starts at that state.
    SciMLBase.set_u!(companion, zero(companion.u))
    state = _StreamingCompanion(companion, copy(prob.u0), false, storage)
    state_ref[] = state
    for _ in integrator
        _advance_streaming_companion!(state, integrator)
    end
    SciMLBase.successful_retcode(integrator.sol) ||
        throw(ErrorException("the forward solve failed with retcode $(integrator.sol.retcode)"))

    # As in InterpolatingMode: the error against P(T⁻), corrected to the saved state.
    endpoint_error = companion.u .- (integrator.u .- state.last_end)
    return _defect_correction_result(
        name, forward, callbacks, integrator.sol, endpoint_error, trajectory;
        verbose = _verbosity(kwargs), t_end = integrator.t
    )
end
