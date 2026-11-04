# Events and semi-explicit DAEs in the adjoint error estimator GlobalAdjoint. An event of the
# forward solve is a jump of P; for the dual it is a jump z ← Sᵀ z across the event, with S the saltation matrix of the event (the linearized map of a
# small error from just before the event to just after it, including the shift of the event
# time). The bridge, the knot terms and the DAE parts of the adjoint solve are in the
# SciMLSensitivity extension.

# The differential and algebraic components of a semi-explicit DAE: a diagonal mass matrix
# with entries 0 and 1 (or the identity for an ODE). A multiple λI of the identity is not
# one: the error identity and the event terms assume M = I.
function _semiexplicit_components(prob, name)
    M = prob.f.mass_matrix
    n = length(prob.u0)
    unsupported = ArgumentError(
        "$name supports only semi-explicit DAEs: a diagonal mass matrix of zeros and ones"
    )
    if M isa LinearAlgebra.UniformScaling
        isone(M.λ) || throw(unsupported)
        return collect(1:n), Int[]
    end
    m = Matrix(M)
    (LinearAlgebra.isdiag(m) && all(x -> iszero(x) || isone(x), LinearAlgebra.diag(m))) ||
        throw(unsupported)
    d = LinearAlgebra.diag(m)
    return findall(isone, d), findall(iszero, d)
end

_rhs(f, iip, u, p, t) = iip ? (du = similar(u); f(du, u, p, t); du) : f(u, p, t)
# The Jacobian of f with ForwardDiff, or with finite differences if ForwardDiff cannot run f
# (the default sensealg also falls back to a numerical VJP then). An error of f itself comes
# again from the finite differences. `verbose` is the verbosity of the solve.
function _rhs_jacobian(f, iip, u, p, t, verbose)
    g = x -> _rhs(f, iip, x, p, t)
    return try
        DifferentiationInterface.jacobian(g, ADTypes.AutoForwardDiff(), u)
    catch err
        err isa InterruptException && rethrow()
        @SciMLMessage(
            "GlobalAdjoint: ForwardDiff cannot differentiate the right-hand side; the " *
                "Jacobian of the DAE at the final time uses finite differences",
            verbose, :global_error_ad_fallback, exception = err, maxlog = 1
        )
        _fd_jacobian(g, u)
    end
end

# The forward events as calls (t, labels). `marks` holds (t, number of events before the
# call), pushed by the forward affect before it records its labels.
function _forward_calls(marks, events)
    return map(eachindex(marks)) do k
        t, start = marks[k]
        stop = k < length(marks) ? marks[k + 1][2] : length(events)
        (t, Tuple(events[(start + 1):stop]))
    end
end

# A parent for the trial runs of the user callbacks (saltation matrix): the options and
# the time axis of the companion integrator, but no side effects (a PresetTimeCallback
# would add a tstop each time).
struct _FrozenParent{O, T, E}
    opts::O
    tdir::T
    dt::T
    tprev::T
    iter::Int
    e::E
end
(parent::_FrozenParent)(t) = parent.e
SciMLBase.add_tstop!(::_FrozenParent, t) = nothing
DiffEqBase.get_tstops_max(parent::_FrozenParent) = typemax(parent.dt)
function _frozen_parent(integrator, e = copy(integrator.u))
    return _FrozenParent(
        integrator.opts, oftype(integrator.t, integrator.tdir), integrator.dt,
        integrator.tprev, integrator.iter, e
    )
end

# Solve the algebraic equations f_a(u) = 0 for the algebraic components (Newton with a
# finite-difference Jacobian; only for the saltation matrix, which needs no high accuracy).
function _project_algebraic(fu, u, alg_idx)
    isempty(alg_idx) && return u
    v = copy(u)
    with_algebraic = x -> (w = copy(v); w[alg_idx] .= x; w)
    for _ in 1:20
        ra = fu(v)[alg_idx]
        LinearAlgebra.norm(ra) <= 1.0e-13 * (1 + LinearAlgebra.norm(v)) && break
        J = _fd_jacobian(x -> fu(with_algebraic(x))[alg_idx], v[alg_idx], ra)
        v[alg_idx] .-= J \ ra
    end
    return v
end

# `g0` is g(x), if the caller has it.
function _fd_jacobian(g, x, g0 = g(x))
    J = similar(x, length(g0), length(x))
    for j in eachindex(x)
        h = sqrt(eps(eltype(x))) * max(1, abs(x[j]))
        w = copy(x)
        w[j] += h
        J[:, j] .= (g(w) .- g0) ./ h
    end
    return J
end

# The user affect of one call, applied to a copy of the true state u and of the parameters.
# Returns the state and the parameters after the affect. `callback` is the trial copy of the
# user callback (`_recorded_forward_solve`): as in SciMLSensitivity, the affect is a
# function of u, p and t, so the trial runs need no copy of the callback at each forward call.
# `args` are the arguments of the forward call after the integrator: for a
# VectorContinuousCallback the vector of the fired components and their signs (without
# them, one built from `labels`).
function _trial_affect(callback, labels, u, t, parent, p, sol, P; args = nothing)
    # Each run gets its own copy of the trial callback, so that the runs of a saltation
    # matrix start from the same callback state (an affect with state in the callback could
    # pass `_check_trial_affect` and still give each finite difference another state). The
    # trial callback never advances, so the copy stays small.
    callback = deepcopy(callback)
    integrator = _TrueStateIntegrator(copy(u), deepcopy(p), t, parent, sol, P)
    _, kind = first(labels)
    if callback isa SciMLBase.ContinuousCallback
        (kind == 1 ? callback.affect! : callback.affect_neg!)(integrator)
    elseif callback isa SciMLBase.VectorContinuousCallback
        callback.affect!(integrator, args === nothing ? _fired_mask(callback, labels) : copy(only(args)))
    else
        callback.affect!(integrator)
    end
    return integrator.u, integrator.p
end

function _fired_mask(callback, labels)
    mask = zeros(Int8, callback.len)
    for (_, component) in labels
        mask[component] = 1
    end
    return mask
end

function _trial_condition(callback, labels, u, t, parent, p, sol, P)
    integrator = _TrueStateIntegrator(copy(u), p, t, parent, sol, P)
    if callback isa SciMLBase.VectorContinuousCallback
        out = similar(u, callback.len)
        callback.condition(out, _true_state_view(u, callback.idxs), t, integrator)
        return out[last(first(labels))]
    end
    return callback.condition(_true_state_view(u, callback.idxs), t, integrator)
end

# The affect of a call and of the calls `tied` to it (each with its trial callback, labels,
# frozen parent and arguments), in order, as a function of the state and the parameters;
# `reinit(u, q)` follows each affect (the reinitialization of a DAE).
function _tied_affects(callback, labels, t, parent, args, tied, sol, P, reinit)
    return function (u, q)
        u, q = _trial_affect(callback, labels, u, t, parent, q, sol, P; args)
        u = reinit(u, q)
        for d in tied
            u, q = _trial_affect(d.callback, d.labels, u, t, d.parent, q, sol, P; d.args)
            u = reinit(u, q)
        end
        return u, q
    end
end

# The saltation matrix of a forward event at t (one affect call), for the differential
# components: S = A + (φ⁺ - A φ⁻) ∇γᵀ / (∇γ ⋅ φ⁻ + γ_t), with A the Jacobian of the affect
# (and the reinitialization of a DAE), γ the event condition and φ the vector field, all
# for the reduced ODE x' = f_d(x, z(x)) of a DAE. Also the offset c = (φ⁺ - A φ⁻) γ₀ /
# (∇γ ⋅ φ⁻ + γ_t) for a condition γ₀ = γ(x_f) that is not zero at the forward event (the
# forward root finding reads an interpolant; for a condition on an algebraic variable that
# interpolant is not consistent): the linearized error after the event is S e⁻ + c.
# A DiscreteCallback does not move its event time: S = A, c = 0. `p` holds the parameters
# before the call. An affect can change them: the algebraic equations after the affect use
# the parameters of the trial run, and φ⁺ uses the state `u_after` (by default P(t⁺);
# another call at the same time comes in between) and the parameters `p_after` after the
# call. A, ∇γ, γ_t and φ⁻ are taken at the state before the call: P(t⁻) for the first call
# at a time, else `u_before`, the state after the calls before it (the calls act in order,
# so φ⁻ of a call is φ⁺ of the call before it). `u_event` is the forward state at the event.
function _event_saltation(
        callback, labels, t, p, parent, P, sol, prob, dif, alg_idx;
        u_event, p_after, u_after = nothing, u_before = nothing, args = nothing, tied = ()
    )
    f = SciMLBase.unwrapped_f(prob.f)
    iip = SciMLBase.isinplace(prob)
    p = deepcopy(p)
    rhs_with = q -> (u -> _rhs(f, iip, u, q, t))
    fu = rhs_with(p)
    Pm = u_before === nothing ? sol(t; continuity = :left) : u_before
    Pp = u_after === nothing ? sol(t; continuity = :right) : u_after
    full = x -> (u = copy(Pm); u[dif] .= x; _project_algebraic(fu, u, alg_idx))
    xm = Pm[dif]
    # The affect of the call and of the `tied` calls after it, each followed by the
    # reinitialization of a DAE.
    affect_all = _tied_affects(
        callback, labels, t, parent, args, tied, sol, P, (u, q) -> _project_algebraic(rhs_with(q), u, alg_idx)
    )
    affected = x -> first(affect_all(full(x), p))[dif]
    A = _fd_jacobian(affected, xm)
    callback isa SciMLBase.DiscreteCallback && return A, zero(xm)
    γ = x -> _trial_condition(callback, labels, full(x), t, parent, p, sol, P)
    ∇γ = vec(_fd_jacobian(x -> [γ(x)], xm))
    h = sqrt(eps(typeof(t))) * max(1, abs(t))
    γt = (_trial_condition(callback, labels, Pm, t + h, parent, p, sol, P) -
        _trial_condition(callback, labels, Pm, t - h, parent, p, sol, P)) / (2h)
    φm = fu(full(xm))[dif]
    φp = rhs_with(p_after)(Pp)[dif]
    jump = (φp .- A * φm) ./ (LinearAlgebra.dot(∇γ, φm) + γt)
    S = A .+ jump * ∇γ'
    # γ₀ at the state where the forward solve evaluated the condition (made consistent): at
    # P(t⁻) it would count again the mismatch m = P(t⁻) - u_f, which has its own term.
    c = jump .* γ(u_event[dif])
    all(isfinite, S) && all(isfinite, c) || return A, zero(xm)
    return S, c
end

# The saltation map as the two products that the dual needs: z ↦ Sᵀz (the dual crosses the
# event) and x ↦ S x (the source of an earlier call at the same time goes through S).
_saltation_maps(S) = (z -> S' * z, x -> S * x)

# The VJP choice of the adjoint method (`_event_vjp`): the two products without S. With
# S = A + j ∇γᵀ and j = (φ⁺ - A φ⁻) / (∇γ ⋅ φ⁻ + γ_t), Sᵀz = Aᵀz + ∇γ (j ⋅ z) and S x = A x + j (∇γ ⋅ x): one
# VJP of the affect for each z (ReverseDiff, which SciMLSensitivity loads; the state is a
# vector of tracked reals, so an affect can write into it), JVPs of the affect (ForwardDiff),
# one gradient of the condition, and γ_t; c = j γ₀. The cost per event does not grow with
# the size of the state, and S is not stored. ODEs only. Returns nothing if the AD fails on
# the user callback (trial products run here, not later in the adjoint solve). The base
# state is the one of `_event_saltation`.
function _saltation_vjp(
        callback, labels, t, p, parent, P, sol, prob;
        u_event, p_after, u_after = nothing, u_before = nothing, args = nothing, tied = ()
    )
    f = SciMLBase.unwrapped_f(prob.f)
    iip = SciMLBase.isinplace(prob)
    p = deepcopy(p)
    Pm = u_before === nothing ? sol(t; continuity = :left) : u_before
    Pp = u_after === nothing ? sol(t; continuity = :right) : u_after
    affect_all = _tied_affects(callback, labels, t, parent, args, tied, sol, P, (u, q) -> u)
    affect = x -> first(affect_all(x, p))
    condition = (x, s) -> _trial_condition(callback, labels, x, s, parent, p, sol, P)
    Atz = z -> DifferentiationInterface.gradient(
        x -> LinearAlgebra.dot(z, affect(collect(x))), ADTypes.AutoReverseDiff(), Pm
    )
    Ax = x -> DifferentiationInterface.derivative(
        h -> affect(Pm .+ h .* x), ADTypes.AutoForwardDiff(), zero(t)
    )
    return try
        Atz(one.(Pm))
        Ax(one.(Pm))
        callback isa SciMLBase.DiscreteCallback && return Atz, Ax, zero(Pm)
        φm = _rhs(f, iip, Pm, p, t)
        φp = _rhs(f, iip, Pp, p_after, t)
        ∇γ = DifferentiationInterface.gradient(
            x -> condition(collect(x), t), ADTypes.AutoReverseDiff(), Pm
        )
        h = sqrt(eps(typeof(t))) * max(1, abs(t))
        γt = (condition(Pm, t + h) - condition(Pm, t - h)) / (2h)
        j = (φp .- Ax(φm)) ./ (LinearAlgebra.dot(∇γ, φm) + γt)
        all(isfinite, j) || return Atz, Ax, zero(Pm)
        (z -> Atz(z) .+ ∇γ .* LinearAlgebra.dot(j, z)),
            (x -> Ax(x) .+ j .* LinearAlgebra.dot(∇γ, x)), j .* condition(u_event, t)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
end

# Merge the windows of several calls at one event time (each window has zero width). They act
# in order, so the bridge is the composition of the saltation maps (`St` applies Sᵀ, `Sx`
# applies S), the first call's e(a⁻) stays, and the source of an earlier call goes through
# the later ones: e(b⁺) = S₂ e₁(b⁺) + e₂(b⁺) on the differential components `dif`.
function _merge_windows(windows, dif)
    merged = Any[]
    for w in sort(windows; by = w -> w.a)
        if !isempty(merged) && w.a == merged[end].b
            m = merged[end]
            eb = copy(w.eb)
            eb[dif] .+= w.Sx(m.eb[dif])
            St1, St2, Sx1, Sx2 = m.St, w.St, m.Sx, w.Sx
            merged[end] = (;
                a = m.a, b = w.b, ea = m.ea, eb,
                St = z -> St1(St2(z)), Sx = x -> Sx2(Sx1(x)),
            )
        else
            push!(merged, w)
        end
    end
    return merged
end

# The forward solve of GlobalAdjoint with its event records: each forward event is a window
# of zero width at its time, with the saltation matrix of the event as the bridge of the dual.
# A record of the state, the parameters and the arguments before each affect serves the
# trial runs, which use one copy of the callbacks for the whole solve (`trial_callbacks`, made
# after the initialization; see `_check_trial_affect`). The affect acts on the integrator state
# u_f, which can differ from the left limit P(t_f⁻) of the dense output (the interpolant of a
# step that the event cuts short); then e⁺ = S (e⁻ + m) + c with m = P(t_f⁻) - u_f and the
# condition offset c (see `_event_saltation`), and the window's pair term z⁺ᵀ (e(b⁺) - S e(a⁻))
# with e(a⁻) = -m and e(b⁺) = c gives the source z⁺ᵀ (S m + c). `alg` is the GlobalAdjoint.
function _recorded_forward_solve(prob, alg, args...; abstol, reltol, kwargs...)
    name = "GlobalAdjoint"
    inner_alg = alg.alg
    event_vjp = _event_vjp(alg)
    supports_events = _supports_events(alg.sensealg)
    verbose = _verbosity(kwargs)
    dif, alg_idx = _semiexplicit_components(prob, name)
    # Without support for events (an adjoint method other than InterpolatingAdjoint), fail
    # before the solve for a DAE, and at the first call that changes the state or the
    # parameters (a call that changes nothing, such as a SavingCallback, is no event).
    supports_events || isempty(alg_idx) || throw(_unsupported_events_error())
    callback_list = _callback_list(_user_callbacks(prob, kwargs))
    events = Tuple{Int, Int}[]
    marks = Tuple{float(eltype(prob.tspan)), Int}[]
    snapshots = Any[]
    # The number of saved states before each call, which places the saved entries at an
    # event time between the calls there (a callback can save before or after its affect).
    saves = Int[]
    # The state after each affect, whether the affect changed the state or the parameters,
    # and then the parameters after it (for `_check_trial_affect`).
    after_states = Any[]
    # Whether an affect at the current time changed the state or the parameters (see
    # `_save_after_events`).
    changed_here = Ref(false)
    after = function (integrator)
        before_call = snapshots[end]
        changed = integrator.u != before_call.u || !_approx_same(before_call.p, integrator.p)
        p = changed ? deepcopy(integrator.p) : nothing
        push!(after_states, (; u = copy(integrator.u), p, changed))
        changed && !supports_events && throw(_unsupported_events_error())
        changed_here[] |= changed
        return nothing
    end
    before = function (integrator, label, args...)
        push!(marks, (integrator.t, length(events)))
        push!(saves, length(integrator.sol.t))
        push!(
            snapshots, (;
                p = deepcopy(integrator.p), parent = _frozen_parent(integrator, zero(integrator.u)),
                u = copy(integrator.u), args = deepcopy(args),
            )
        )
        return after
    end
    forward_prob, forward_callbacks = _forward_setup(prob, callback_list, events, before)
    # `init` and `solve!`, so that the parameters at the end are available.
    integrator = SciMLBase.init(
        forward_prob, inner_alg, args...;
        abstol, reltol,
        callback = SciMLBase.CallbackSet(forward_callbacks, _save_after_events(changed_here)),
        merge(_without_callback(kwargs), _DENSE_SOLVE_KWARGS)...
    )
    # A change by an initialization has its own term (e0 below).
    changed_here[] = false
    # The state after the initialization of the callbacks: the error at t₀ is the difference
    # to the start of the dense output (zero unless an initialization changes the state and
    # the solution does not save the result).
    u_init = copy(integrator.u)
    trial_callbacks = deepcopy(callback_list)
    sol = SciMLBase.solve!(integrator)
    _validate_estimation_solution(sol, name)
    knots, jumps = _knot_jumps(sol)
    P = _IntervalInterp(sol, [knots; sol.t[end]], 1)
    e0 = u_init .- sol(sol.t[1]; continuity = :right)
    calls = _forward_calls(marks, events)
    params = _parameter_record(calls, snapshots, deepcopy(integrator.p))
    zero_e = zero(prob.u0)
    # Several calls at one time act in order: the state and the parameters after a call are
    # those before the next one; only the first call has the cut-step mismatch. After the
    # last call, the dense output starts at P(t⁺), the state after the calls
    # (`_save_after_events`). An affect that changes neither the state nor the parameters (a
    # SavingCallback) has S = I and no source, to first order, and needs no trial runs.
    # A DiscreteCallback call right after a continuous call at the same time acts at the step
    # end that the event makes, in every solve, so it moves with the event time: the window
    # of the continuous call takes the affects of those calls too (one saltation matrix, with
    # the vector field after all of them), and their own windows are the identity.
    discrete(k) = _is_discrete(first(last(calls[k])))
    tied_to = zeros(Int, length(calls))
    last_tied = collect(eachindex(calls))
    for k in eachindex(calls)
        discrete(k) && continue
        j = k
        while j < length(calls) && first(calls[j + 1]) == first(calls[k]) && discrete(j + 1)
            j += 1
            tied_to[j] = k
        end
        last_tied[k] = j
    end
    windows = map(enumerate(calls)) do (k, (tf, labels))
        (; p, parent, args) = snapshots[k]
        callback = trial_callbacks[first(first(labels))]
        u_f = snapshots[k].u
        # The error that the forward affect makes (zero for an affect that changes nothing).
        affect_error = after_states[k].changed ?
            _check_trial_affect(callback, labels, tf, snapshots[k], after_states[k], P, sol) :
            zero_e
        if tied_to[k] != 0
            # (Its affect error is added here, not passed through the affects after it in the
            # group: a term of the size of Δ times (A - I).) Not a variable `eb`: `bridge`
            # below assigns that name, so it would capture this one in a box.
            source = copy(zero_e)
            source[dif] .-= affect_error[dif]
            return (;
                a = tf, b = tf, ea = copy(zero_e), eb = source, St = identity,
                Sx = identity, saves = saves[k], group_saves = [saves[k]], partial = nothing,
                after_states[k].changed,
            )
        end
        group = k:last_tied[k]
        first_at_time = k == 1 || first(calls[k - 1]) != tf
        u_before = first_at_time ? nothing : u_f
        # The bridge of the first `m` calls of the group: the saltation of their affects, with
        # the state and the parameters after them (those of the next call, or P(t⁺) and the
        # parameters after the event time), and the source term eb.
        function bridge(m)
            part = group[1:m]
            # (Only the calls that changed something: the others have A = I, as everywhere.
            # The affect of a SavingCallback cannot take the dual numbers of AD.)
            tied = [
                (;
                    callback = trial_callbacks[first(first(calls[i][2]))], labels = calls[i][2],
                    snapshots[i].parent, snapshots[i].args,
                ) for i in part[2:end] if after_states[i].changed
            ]
            changed = any(i -> after_states[i].changed, part)
            next = last(part) < length(calls) && first(calls[last(part) + 1]) == tf ?
                snapshots[last(part) + 1] : nothing
            p_after = next === nothing ? _parameters_at(params, tf, :right) : next.p
            u_after = next === nothing ? nothing : next.u
            vjp = event_vjp && isempty(alg_idx) && changed
            maps = if !changed
                identity, identity, zero(u_f[dif])
            elseif vjp
                _saltation_vjp(
                    callback, labels, tf, p, parent, P, sol, prob;
                    u_event = u_f, p_after, u_after, u_before, args, tied
                )
            else
                nothing
            end
            if maps === nothing
                vjp && @SciMLMessage(
                    "GlobalAdjoint uses the finite-difference saltation matrix at the event at " *
                        "t = $tf: ReverseDiff and ForwardDiff cannot differentiate the callback. " *
                        "Set `autojacvec = false` in the `sensealg` to use finite differences " *
                        "without this warning.",
                    verbose, :global_error_ad_fallback, maxlog = 1
                )
                S, c = _event_saltation(
                    callback, labels, tf, p, parent, P, sol, prob, dif, alg_idx;
                    u_event = u_f, p_after, u_after, u_before, args, tied
                )
                St, Sx = _saltation_maps(S)
            else
                St, Sx, c = maps
            end
            eb = copy(zero_e)
            eb[dif] .= c .- affect_error[dif]
            return (; St, Sx, eb)
        end
        # (Not `(; St, Sx, eb) = ...`: `bridge` assigns those names, so it would write to them.)
        whole = bridge(length(group))
        # A saved entry between the calls of a group (a callback that saves before or after
        # its affect) is the state after the first m calls: its own bridge (`_entry_bridge`).
        # Only an entry that an estimate ends at needs it, so it is computed then.
        group_saves = saves[group]
        partial = if length(group) == 1
            nothing
        else
            partials = Dict{Int, Any}()
            m -> get!(() -> bridge(m), partials, m)
        end
        ea = first_at_time ? u_f .- sol(tf; continuity = :left) : copy(zero_e)
        (;
            a = tf, b = tf, ea, whole.eb, whole.St, whole.Sx, saves = saves[k], group_saves,
            partial, after_states[k].changed,
        )
    end
    return (;
        prob, sol, knots, jumps, windows = _merge_windows(windows, dif), calls = windows, dif,
        alg_idx, params, e0, verbose,
    )
end

# The trial runs assume that the affect is a function of the state, the parameters and the
# time. Check it at each call that changes them: the affect, run again on the trial copy at
# the state `before` the call, must give the state and the parameters `after` it in the
# forward solve. An affect that reads the derivative or the interpolant gets values in the
# trial run (f(u) and the stored solution) that differ from the forward ones by about the
# local error, so the tolerance is 10 times the forward tolerances; an affect that reads a
# counter in the callback changes the result by much more and fails here. Returns that
# difference Δ = u⁺ - affect(u_f) of the state: an error that the forward affect makes, so
# e⁺ = S e⁻ + c - Δ (the window's source term).
function _check_trial_affect(callback, labels, t, before, after, P, sol)
    u, p = _trial_affect(
        callback, labels, before.u, t, before.parent, before.p, sol, P; args = before.args
    )
    opts = before.parent.opts
    rtol = 10 * maximum(opts.reltol)
    atol = 10 * maximum(opts.abstol)
    (_approx_same(u, after.u; rtol, atol) && _approx_same(p, after.p; rtol, atol)) || throw(
        ArgumentError(
            "GlobalAdjoint: the affect of the event at t = $t gives another result when it " *
                "runs again at the same state, parameters and time. GlobalAdjoint needs an " *
                "affect that is a function of these, as SciMLSensitivity's adjoints do, not " *
                "of state kept in the callback (such as a counter). GlobalErrorEstimation " *
                "supports such affects."
        )
    )
    return after.u .- u
end

# Whether a and b are equal up to the tolerances (exactly equal by default): numbers and
# arrays of numbers by `isapprox`, other values by `isequal`, or else element by element (an
# array) or field by field (a struct without `==`, for which `isequal` is false for a copy,
# as for the parameters before and after an affect). An error in the comparison counts as a
# difference.
function _approx_same(a, b; rtol = 0, atol = 0)
    return try
        _approx_same_values(a, b, rtol, atol)
    catch err
        err isa InterruptException && rethrow()
        false
    end
end

function _approx_same_values(a, b, rtol, atol)
    isequal(a, b) && return true
    numeric = x -> x isa Union{Number, AbstractArray{<:Number}}
    numeric(a) && numeric(b) && return isapprox(a, b; rtol, atol, nans = true)
    a isa AbstractArray && b isa AbstractArray && return size(a) == size(b) &&
        all(_approx_same_values(x, y, rtol, atol) for (x, y) in zip(a, b))
    T = typeof(a)
    return T === typeof(b) && isstructtype(T) && nfields(a) > 0 &&
        all(i -> _approx_same_values(getfield(a, i), getfield(b, i), rtol, atol), 1:nfields(a))
end

function _unsupported_events_error()
    return ArgumentError(
        "GlobalAdjoint supports DAEs and events that change the state or the parameters only " *
            "with an InterpolatingAdjoint sensealg (the default)"
    )
end

# The dense output of the step after an event starts from the last saved state. After the
# calls at an event time that changed the state or the parameters, this callback saves the
# state, if no callback saved it, so that P(t⁺) is the state after the calls. (Else P on
# that step starts from the state before the affect, far from the true solution, and the
# identity is a cancellation of large terms.) It comes after the user callbacks, and so
# after the reinitialization of a DAE.
function _save_after_events(changed_here)
    return SciMLBase.DiscreteCallback(
        (u, t, integrator) -> changed_here[],
        function (integrator)
            changed_here[] = false
            sol = integrator.sol
            sol.t[end] == integrator.t && sol.u[end] == integrator.u ||
                SciMLBase.savevalues!(integrator, true)
            SciMLBase.u_modified!(integrator, false)
            return nothing
        end;
        save_positions = (false, false), initializealg = SciMLBase.NoInit()
    )
end

# The calls at the time of the saved entry `index` that act before it was saved, in order.
function _calls_before_entry(data, index)
    t = data.sol.t[index]
    return [w for w in data.calls if w.b == t && w.saves < index]
end

# The bridge (St, Sx, eb) of the call window `w` for the saved entry `index`: the window's own,
# or, for an entry between the calls of its group, the bridge of the calls before the entry.
function _entry_bridge(w, index)
    m = count(<(index), w.group_saves)
    return m < length(w.group_saves) ? w.partial(m) : w
end

# The parameters of the forward solve on each interval between its event times (an affect
# can change them). The value before each call is in its snapshot; the value after the last
# event is the one of the integrator at the end. A change that a DAE reinitialization makes
# after an affect is also included, because the next snapshot comes after it.
function _parameter_record(calls, snapshots, p_end)
    first_call = [k for k in eachindex(calls) if k == 1 || first(calls[k]) != first(calls[k - 1])]
    times = [first(calls[k]) for k in first_call]
    values = [[snapshots[k].p for k in first_call]; [p_end]]
    return (; times, values = map(identity, values))
end

# The parameters at t, on the left or right side of an event time.
function _parameters_at(params, t, side)
    (; times, values) = params
    k = side === :left ? searchsortedfirst(times, t) : searchsortedlast(times, t) + 1
    return values[k]
end

# The same data for a solution without callbacks (no events), such as the fixed-grid solve of
# StepGridRefinement.
function _unrecorded_data(prob, sol, verbose)
    dif, alg_idx = _semiexplicit_components(prob, "GlobalAdjoint")
    knots, jumps = _knot_jumps(sol)
    params = (; times = eltype(sol.t)[], values = [sol.prob.p])
    return (;
        prob, sol, knots, jumps, windows = Any[], calls = Any[], dif, alg_idx, params,
        e0 = zero(sol.u[1]), verbose,
    )
end
