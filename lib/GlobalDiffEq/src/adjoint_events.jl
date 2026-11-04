# Events and semi-explicit DAEs in GlobalAdjoint. At a forward event the dual jumps by
# z ← Sᵀz, with S the saltation matrix: the linearized map of an error from just before the
# event to just after it, with the shift of the event time. The bridge, the knot terms and
# the DAE parts of the adjoint solve are in the SciMLSensitivity extension.

# The differential and algebraic indices of a semi-explicit DAE (a diagonal mass matrix of
# zeros and ones). Reject λI with λ ≠ 1: the error identity and the event terms assume
# M = I.
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
# Fall back to finite differences if ForwardDiff cannot run f, as the default sensealg falls
# back to a numerical VJP. If f itself throws, `_fd_jacobian` throws the error again.
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

# The forward events grouped into calls (t, labels). `marks[k]` is (t, number of events
# before call k).
function _forward_calls(marks, events)
    return map(eachindex(marks)) do k
        t, start = marks[k]
        stop = k < length(marks) ? marks[k + 1][2] : length(events)
        (t, Tuple(events[(start + 1):stop]))
    end
end

# The parent integrator of the trial runs: the options and time axis of the companion
# integrator, but no side effects (a PresetTimeCallback would add a tstop at each run).
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

# Newton for f_a(u) = 0 with a finite-difference Jacobian. Only the saltation matrix uses
# it, and it needs no high accuracy.
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

# Run the affect of one call on copies of u and p. As in SciMLSensitivity, the affect is a
# function of u, p and t, so one trial copy of the callback serves all calls. For a
# VectorContinuousCallback, `args` holds the fired components and their signs; without it,
# `_fired_mask` builds them from `labels`.
function _trial_affect(callback, labels, u, t, parent, p, sol, P; args = nothing)
    # A new copy for each run, so that all finite differences of one saltation matrix start
    # from the same callback state. The trial callback never advances, so the copy is small.
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

# The affects of a call and of its `tied` calls, in order, as a map (u, q) ↦ (u⁺, q⁺).
# `reinit(u, q)` runs after each affect (the reinitialization of a DAE).
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

# The saltation matrix of a forward event on the differential components,
#     S = A + (φ⁺ - A φ⁻) ∇γᵀ / (∇γ ⋅ φ⁻ + γ_t),
# with A the Jacobian of the affects (and the DAE reinitialization), γ the event condition
# and φ the vector field of the reduced ODE x' = f_d(x, z(x)). Also the offset
# c = (φ⁺ - A φ⁻) γ₀ / (∇γ ⋅ φ⁻ + γ_t), so that e⁺ = S e⁻ + c. Here γ₀ = γ(`u_event`) can be
# nonzero: the forward root finding reads an interpolant, which is not consistent for a
# condition on an algebraic variable. A DiscreteCallback does not move its time, so S = A
# and c = 0.
# A, ∇γ, γ_t and φ⁻ use the state before the call: P(t⁻), or `u_before` for a later call at
# the same time. φ⁺ uses the state `u_after` (default P(t⁺)) and the parameters `p_after`.
function _event_saltation(
        callback, labels, t, p, parent, P, sol, prob, dif, alg_idx;
        u_event, p_after, u_after = nothing, u_before = nothing, args = nothing, tied = (),
        verbose = DEVerbosity()
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
    # γ₀ at the forward event state, not at P(t⁻): the mismatch m = P(t⁻) - u_f has its own
    # source term.
    c = jump .* γ(u_event[dif])
    all(isfinite, S) && all(isfinite, c) || return (_grazing_message(t, verbose); (A, zero(xm)))
    return S, c
end

# At a grazing event ∇γ ⋅ φ⁻ + γ_t is zero, so the time-shift term is not finite and S
# falls back to A.
function _grazing_message(t, verbose)
    @SciMLMessage(
        "GlobalAdjoint: the event at t = $t is grazing (the vector field is tangent to the " *
            "event surface), so its saltation matrix has no time-shift term and the estimate " *
            "is less accurate there.",
        verbose, :global_error_grazing_event
    )
    return nothing
end

# The two products with S that the dual needs: z ↦ Sᵀz for the jump of the dual, and
# x ↦ S x to carry the source of an earlier call at the same time through this call.
_saltation_maps(S) = (z -> S' * z, x -> S * x)

# The matrix-free form of `_event_saltation` (see `_event_vjp`), for ODEs only. With
# j = (φ⁺ - A φ⁻) / (∇γ ⋅ φ⁻ + γ_t):
#     Sᵀz = Aᵀz + ∇γ (j ⋅ z),   S x = A x + j (∇γ ⋅ x),   c = j γ₀.
# Aᵀz uses ReverseDiff (a tracked-real state lets an affect write into it); A x uses
# ForwardDiff. The cost per event does not grow with the state size. Returns nothing if the
# AD fails on the user callback: the trial products below catch that here, not in the
# adjoint solve.
function _saltation_vjp(
        callback, labels, t, p, parent, P, sol, prob;
        u_event, p_after, u_after = nothing, u_before = nothing, args = nothing, tied = (),
        verbose = DEVerbosity()
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
        all(isfinite, j) || return (_grazing_message(t, verbose); (Atz, Ax, zero(Pm)))
        (z -> Atz(z) .+ ∇γ .* LinearAlgebra.dot(j, z)),
            (x -> Ax(x) .+ j .* LinearAlgebra.dot(∇γ, x)), j .* condition(u_event, t)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
end

# Merge the zero-width windows of the calls at one event time. The calls act in order, so
# the bridge composes their saltation maps, the first call's e(a⁻) stays, and an earlier
# source goes through the later maps: e(b⁺) = S₂ e₁(b⁺) + e₂(b⁺) on `dif`.
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

# The forward solve of GlobalAdjoint with its event records. Each call is a zero-width
# window whose bridge is its saltation matrix. The trial runs use the snapshots before each
# affect and one copy of the callbacks, made after the initialization
# (`_check_trial_affect`). The affect acts on the integrator state u_f, which can differ
# from P(t_f⁻) if the event cuts the step short. With m = P(t_f⁻) - u_f and the offset c of
# `_event_saltation`,
#     e⁺ = S (e⁻ + m) + c,
# so the pair term z⁺ᵀ (e(b⁺) - S e(a⁻)) with e(a⁻) = -m and e(b⁺) = c gives z⁺ᵀ (S m + c).
function _recorded_forward_solve(prob, alg, args...; abstol, reltol, kwargs...)
    name = "GlobalAdjoint"
    inner_alg = alg.alg
    event_vjp = _event_vjp(alg)
    supports_events = _supports_events(alg.sensealg)
    verbose = _verbosity(kwargs)
    dif, alg_idx = _semiexplicit_components(prob, name)
    # Without event support (a sensealg other than InterpolatingAdjoint), fail before the
    # solve for a DAE, and at the first call that changes u or p. A call that changes
    # nothing (such as a SavingCallback) is not an event.
    supports_events || isempty(alg_idx) || throw(_unsupported_events_error())
    callback_list = _callback_list(_user_callbacks(prob, kwargs))
    events = Tuple{Int, Int}[]
    marks = Tuple{float(eltype(prob.tspan)), Int}[]
    snapshots = Any[]
    # The saved-state count before each call places the saved entries at an event time
    # between the calls at that time (a callback can save before or after its affect).
    saves = Int[]
    # For each affect: the state after it, whether it changed u or p, and p after it if it
    # did (for `_check_trial_affect`).
    after_states = Any[]
    # Whether an affect at the current time changed u or p (for `_save_after_events`).
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
    # `init` and `solve!`, so that the integrator holds the final parameters.
    integrator = SciMLBase.init(
        forward_prob, inner_alg, args...;
        abstol, reltol,
        callback = SciMLBase.CallbackSet(forward_callbacks, _save_after_events(changed_here)),
        merge(_without_callback(kwargs), _DENSE_SOLVE_KWARGS)...
    )
    # A change by an initialization goes into e0, not into an event window.
    changed_here[] = false
    # e0 is the state after the callback initialization minus the start of the dense
    # output. It is zero unless an initialization changes the state without a save.
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
    # Calls at one time act in order: the state and parameters after a call are those before
    # the next one, and only the first call has the cut-step mismatch. After the last call,
    # the dense output starts at P(t⁺) (`_save_after_events`). A call that changes neither u
    # nor p has S = I and no source to first order, so it needs no trial runs.
    # A DiscreteCallback call right after a continuous call at the same time acts at the
    # step end that the event makes in any nearby solve, so it moves with the event time.
    # Hence the continuous call's window takes the affects of these tied calls too (one
    # saltation matrix, with φ⁺ after all of them), and their own windows are the identity.
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
        # Δ of `_check_trial_affect`, or zero for an affect that changes nothing.
        affect_error = after_states[k].changed ?
            _check_trial_affect(callback, labels, tf, snapshots[k], after_states[k], P, sol) :
            zero_e
        if tied_to[k] != 0
            # Add Δ after the group, not through its later affects; this drops a term of
            # size Δ (A - I). A new name, so that the closure `bridge` does not capture this
            # variable in a box.
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
        # The bridge of the first `m` calls of the group. The state and parameters after
        # them are those of the next call, or P(t⁺) and the parameters after the event time.
        function bridge(m)
            part = group[1:m]
            # Only the calls that changed something: the others have A = I, and the affect
            # of a SavingCallback cannot take AD dual numbers.
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
                    u_event = u_f, p_after, u_after, u_before, args, tied, verbose
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
                    u_event = u_f, p_after, u_after, u_before, args, tied, verbose
                )
                St, Sx = _saltation_maps(S)
            else
                St, Sx, c = maps
            end
            eb = copy(zero_e)
            eb[dif] .= c .- affect_error[dif]
            return (; St, Sx, eb)
        end
        # Not `(; St, Sx, eb) = ...`: the closure `bridge` would then write to them.
        whole = bridge(length(group))
        # A saved entry between the calls of a group needs the bridge of the calls before
        # it (`_entry_bridge`). Only an entry where an estimate ends uses it, so compute it
        # lazily.
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

# The trial runs assume that the affect is a function of u, p and t. Check this at each
# call that changes u or p: a rerun at the state `before` the call must give the state
# `after` it. An affect that reads the derivative or the interpolant differs by about the
# local error, so the tolerance is 10 times the solver tolerances. An affect that reads a
# counter fails. Returns the forward affect error Δ = u⁺ - affect(u_f), so that
# e⁺ = S e⁻ + c - Δ.
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

# Equality up to the tolerances (exact by default). Compare element by element or field by
# field if `isequal` fails: a struct without `==` is not `isequal` to its copy. A comparison
# that throws counts as a difference.
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

# The dense output of the step after an event starts from the last saved state. Save the
# state after the calls that changed u or p, so that P(t⁺) is that state. Otherwise P starts
# from the state before the affect, and the error identity cancels large terms. This
# callback comes after the user callbacks, so after the DAE reinitialization.
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

# The calls at the time of saved entry `index` that act before the save, in order.
function _calls_before_entry(data, index)
    t = data.sol.t[index]
    return [w for w in data.calls if w.b == t && w.saves < index]
end

# The bridge (St, Sx, eb) of call window `w` for saved entry `index`. An entry between the
# calls of a group gets the bridge of the calls before it.
function _entry_bridge(w, index)
    m = count(<(index), w.group_saves)
    return m < length(w.group_saves) ? w.partial(m) : w
end

# The parameters on each interval between event times, since an affect can change them.
# The next snapshot comes after a DAE reinitialization, so it includes any change from it.
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

# The same data for a solution without events, such as the fixed-grid solve of
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
