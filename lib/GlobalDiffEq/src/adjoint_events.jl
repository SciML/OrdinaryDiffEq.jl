# Events and semi-explicit DAEs in the adjoint error estimators (GlobalAdjoint, and
# CompanionAdjoint). An event of the forward solve is a jump of P; for the dual it is a jump
# z ← Sᵀ z across the event, with S the saltation matrix of the event (the linearized map of a
# small error from just before the event to just after it, including the shift of the event
# time). The bridge, the knot terms and the DAE parts of the adjoint solve are in the
# SciMLSensitivity extension.

# The differential and algebraic components of a semi-explicit DAE: a diagonal mass matrix
# with entries 0 and 1 (or the identity for an ODE).
function _semiexplicit_components(prob, name)
    M = prob.f.mass_matrix
    n = length(prob.u0)
    M isa LinearAlgebra.UniformScaling && return collect(1:n), Int[]
    m = Matrix(M)
    (LinearAlgebra.isdiag(m) && all(x -> iszero(x) || isone(x), LinearAlgebra.diag(m))) || throw(
        ArgumentError(
            "$name supports only semi-explicit DAEs: a diagonal mass matrix of zeros and ones"
        )
    )
    d = LinearAlgebra.diag(m)
    return findall(isone, d), findall(iszero, d)
end

_rhs(f, iip, u, p, t) = iip ? (du = similar(u); f(du, u, p, t); du) : f(u, p, t)
_rhs_jacobian(f, iip, u, p, t) = DifferentiationInterface.jacobian(
    x -> _rhs(f, iip, x, p, t), ADTypes.AutoForwardDiff(), u
)

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
    u = copy(u)
    for _ in 1:20
        ra = fu(u)[alg_idx]
        LinearAlgebra.norm(ra) <= 1.0e-13 * (1 + LinearAlgebra.norm(u)) && break
        J = similar(u, length(alg_idx), length(alg_idx))
        for (j, k) in enumerate(alg_idx)
            h = sqrt(eps(eltype(u))) * max(1, abs(u[k]))
            w = copy(u)
            w[k] += h
            J[:, j] .= (fu(w)[alg_idx] .- ra) ./ h
        end
        u[alg_idx] .-= J \ ra
    end
    return u
end

function _fd_jacobian(g, x)
    g0 = g(x)
    J = similar(x, length(g0), length(x))
    for j in eachindex(x)
        h = sqrt(eps(eltype(x))) * max(1, abs(x[j]))
        w = copy(x)
        w[j] += h
        J[:, j] .= (g(w) .- g0) ./ h
    end
    return J
end

# The user affect of one call, applied to a copy of the true state u (and of the callback,
# which can hold state, such as the next time of a PresetTimeCallback, and of the
# parameters). Returns the state and the parameters after the affect.
function _trial_affect(callback, labels, u, t, parent, p, sol, P)
    callback = deepcopy(callback)
    integrator = _TrueStateIntegrator(copy(u), deepcopy(p), t, parent, sol, P)
    _, kind = first(labels)
    if callback isa SciMLBase.ContinuousCallback
        (kind == 1 ? callback.affect! : callback.affect_neg!)(integrator)
    elseif callback isa SciMLBase.VectorContinuousCallback
        callback.affect!(integrator, kind)
    else
        callback.affect!(integrator)
    end
    return integrator.u, integrator.p
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

# The saltation matrix of a forward event at t (one affect call), for the differential
# components: S = A + (φ⁺ - A φ⁻) ∇γᵀ / (∇γ ⋅ φ⁻ + γ_t), with A the Jacobian of the affect
# (and the reinitialization of a DAE), γ the event condition and φ the vector field, all
# for the reduced ODE x' = f_d(x, z(x)) of a DAE. Also the offset c = (φ⁺ - A φ⁻) γ₀ /
# (∇γ ⋅ φ⁻ + γ_t) for a condition γ₀ = γ(x_f) that is not zero at the forward event (the
# forward root finding reads an interpolant; for a condition on an algebraic variable that
# interpolant is not consistent): the linearized error after the event is S e⁻ + c.
# A DiscreteCallback does not move its event time: S = A, c = 0. For CompanionAdjoint any S
# keeps the estimate exact (a poor S only makes the terms larger); GlobalAdjoint needs S and c.
# `p` holds the parameters before the call. An affect can change them: the algebraic
# equations after the affect use the parameters of the trial run, and φ⁺ uses the state
# `u_after` and the parameters `p_after` after the call (by default P(t⁺) and the parameters
# of the trial run; another call at the same time comes in between).
function _event_saltation(
        callback, labels, t, p, parent, P, sol, prob, dif, alg_idx;
        u_event = nothing, p_after = nothing, u_after = nothing
    )
    f = SciMLBase.unwrapped_f(prob.f)
    iip = SciMLBase.isinplace(prob)
    p = deepcopy(p)
    rhs_with = q -> (u -> _rhs(f, iip, u, q, t))
    fu = rhs_with(p)
    Pm = sol(t; continuity = :left)
    Pp = u_after === nothing ? sol(t; continuity = :right) : u_after
    full = x -> (u = copy(Pm); u[dif] .= x; _project_algebraic(fu, u, alg_idx))
    xm = Pm[dif]
    affected = function (x)
        u, q = _trial_affect(callback, labels, full(x), t, parent, p, sol, P)
        return _project_algebraic(rhs_with(q), u, alg_idx)[dif]
    end
    A = _fd_jacobian(affected, xm)
    callback isa SciMLBase.DiscreteCallback && return A, zero(xm)
    γ = x -> _trial_condition(callback, labels, full(x), t, parent, p, sol, P)
    ∇γ = vec(_fd_jacobian(x -> [γ(x)], xm))
    h = sqrt(eps(typeof(t))) * max(1, abs(t))
    γt = (_trial_condition(callback, labels, Pm, t + h, parent, p, sol, P) -
        _trial_condition(callback, labels, Pm, t - h, parent, p, sol, P)) / (2h)
    φm = fu(full(xm))[dif]
    p_after = p_after === nothing ?
        last(_trial_affect(callback, labels, full(xm), t, parent, p, sol, P)) : p_after
    φp = rhs_with(p_after)(Pp)[dif]
    jump = (φp .- A * φm) ./ (LinearAlgebra.dot(∇γ, φm) + γt)
    S = A .+ jump * ∇γ'
    # γ₀ at the state where the forward solve evaluated the condition (made consistent): at
    # P(t⁻) it would count again the mismatch m = P(t⁻) - u_f, which has its own term.
    c = jump .* γ(u_event === nothing ? xm : u_event[dif])
    all(isfinite, S) && all(isfinite, c) || return A, zero(xm)
    return S, c
end

# The saltation map as the two products that the dual needs: z ↦ Sᵀz (the dual crosses the
# event) and x ↦ S x (the source of an earlier call at the same time goes through S).
_saltation_maps(S) = (z -> S' * z, x -> S * x)

# The VJP choice of the adjoint method (`_event_vjp`): the two products without S. With S = A + j ∇γᵀ and
# j = (φ⁺ - A φ⁻) / (∇γ ⋅ φ⁻ + γ_t), Sᵀz = Aᵀz + ∇γ (j ⋅ z) and S x = A x + j (∇γ ⋅ x): one
# VJP of the affect for each z (ReverseDiff, which SciMLSensitivity loads; the state is a
# vector of tracked reals, so an affect can write into it), JVPs of the affect (ForwardDiff),
# one gradient of the condition, and γ_t; c = j γ₀. The cost per event does not grow with
# the size of the state, and S is not stored. ODEs only. Returns nothing if the AD fails on
# the user callback (a trial product runs here, not later in the adjoint solve).
function _saltation_vjp(
        callback, labels, t, p, parent, P, sol, prob;
        u_event, p_after = nothing, u_after = nothing
    )
    f = SciMLBase.unwrapped_f(prob.f)
    iip = SciMLBase.isinplace(prob)
    p = deepcopy(p)
    Pm = sol(t; continuity = :left)
    Pp = u_after === nothing ? sol(t; continuity = :right) : u_after
    affect = x -> first(_trial_affect(callback, labels, x, t, parent, p, sol, P))
    condition = (x, s) -> _trial_condition(callback, labels, x, s, parent, p, sol, P)
    Atz = z -> DifferentiationInterface.gradient(
        x -> LinearAlgebra.dot(z, affect(collect(x))), ADTypes.AutoReverseDiff(), Pm
    )
    Ax = x -> DifferentiationInterface.derivative(
        h -> affect(Pm .+ h .* x), ADTypes.AutoForwardDiff(), zero(t)
    )
    return try
        Atz(one.(Pm))
        callback isa SciMLBase.DiscreteCallback && return Atz, Ax, zero(Pm)
        φm = _rhs(f, iip, Pm, p, t)
        p_after = p_after === nothing ? last(_trial_affect(callback, labels, Pm, t, parent, p, sol, P)) : p_after
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

# Merge windows that overlap (events closer together than their windows, as in a bouncing
# ball, or several events at one time): the merged window spans both, starts with the
# earlier e(a⁻), ends with the later e(b⁺), and its bridge is the composition of the two
# saltation maps, in time order (`St` applies Sᵀ, `Sx` applies S). Two windows of zero width at one time (GlobalAdjoint,
# several calls at one event time) act in order, so the source of the first one goes through
# the second: e(b⁺) = S₂ e₁(b⁺) + e₂(b⁺) on the differential components `dif`.
function _merge_windows(windows; dif = nothing)
    merged = Any[]
    for w in sort(windows; by = w -> w.a)
        if !isempty(merged) && w.a <= merged[end].b
            m = merged[end]
            eb = if w.a == w.b == m.b
                d = something(dif, eachindex(w.eb))
                e = copy(w.eb)
                e[d] .+= w.Sx(m.eb[d])
                e
            else
                w.b > m.b ? w.eb : m.eb
            end
            St1, St2, Sx1, Sx2 = m.St, w.St, m.Sx, w.Sx
            merged[end] = (;
                a = m.a, b = max(m.b, w.b), ea = m.ea, eb,
                St = z -> St1(St2(z)), Sx = x -> Sx2(Sx1(x)), paired = m.paired || w.paired,
            )
        else
            push!(merged, w)
        end
    end
    return merged
end

# The forward solve of GlobalAdjoint with its event records: each forward event is a window
# of zero width at its time, with the saltation matrix of the event as the bridge of the dual
# (the linearized limit of the CompanionAdjoint bridge). A snapshot of the integrator before
# each affect serves the trial runs of the callbacks. The affect acts on the integrator state
# u_f, which can differ from the left limit P(t_f⁻) of the dense output (the interpolant of a
# step that the event cuts short); then e⁺ = S (e⁻ + m) + c with m = P(t_f⁻) - u_f and the
# condition offset c (see `_event_saltation`), and the window's pair term z⁺ᵀ (e(b⁺) - S e(a⁻))
# with e(a⁻) = -m and e(b⁺) = c gives the source z⁺ᵀ (S m + c).
function _recorded_forward_solve(
        prob, inner_alg, name, args...; abstol, reltol, event_vjp = false, kwargs...
    )
    dif, alg_idx = _semiexplicit_components(prob, name)
    callback_list = _callback_list(_user_callbacks(prob, kwargs))
    trial_list = deepcopy(callback_list)
    events = Tuple{Int, Int}[]
    marks = Tuple{eltype(prob.tspan), Int}[]
    snapshots = Any[]
    # The number of saved states before each call, which places the saved entries at an
    # event time between the calls there (a callback can save before or after its affect).
    saves = Int[]
    # The state after each affect: the solution does not save it when the callback saves only
    # before its affect, and then the dense output after the event starts from the state
    # before it.
    after_states = Any[]
    after = integrator -> (push!(after_states, copy(integrator.u)); nothing)
    before = function (integrator)
        push!(marks, (integrator.t, length(events)))
        push!(saves, length(integrator.sol.t))
        push!(snapshots, (deepcopy(integrator.p), _frozen_parent(integrator, zero(integrator.u)), copy(integrator.u)))
        return after
    end
    forward_prob, forward_callbacks = _forward_setup(prob, callback_list, events, before)
    # `init` and `solve!`, so that the parameters at the end are available.
    integrator = SciMLBase.init(
        forward_prob, inner_alg, args...;
        abstol, reltol, callback = forward_callbacks,
        merge(_without_callback(kwargs), _DENSE_SOLVE_KWARGS)...
    )
    sol = SciMLBase.solve!(integrator)
    _validate_estimation_solution(sol, name)
    knots, jumps = _knot_jumps(sol)
    P = _IntervalInterp(sol, [knots; sol.t[end]], 1)
    calls = _forward_calls(marks, events)
    params = _parameter_record(calls, snapshots, deepcopy(integrator.p))
    zero_e = zero(prob.u0)
    # Several calls at one time act in order: the state and the parameters after a call are
    # those before the next one; only the first call has the cut-step mismatch. After the
    # last call, the dense output starts at P(t⁺), which is not the state after the affect
    # when the callback does not save that state: e(t⁺) has the offset u⁺ - P(t⁺) too.
    windows = map(enumerate(calls)) do (k, (tf, labels))
        p, parent, u_f = snapshots[k]
        next = k < length(calls) && first(calls[k + 1]) == tf ? snapshots[k + 1] : nothing
        first_at_time = k == 1 || first(calls[k - 1]) != tf
        callback = trial_list[first(first(labels))]
        p_after = next === nothing ? _parameters_at(params, tf, :right) : first(next)
        u_after = next === nothing ? nothing : last(next)
        vjp = event_vjp && isempty(alg_idx)
        maps = vjp ? _saltation_vjp(
                callback, labels, tf, p, parent, P, sol, prob; u_event = u_f, p_after, u_after
            ) : nothing
        if maps === nothing
            vjp && @warn(
                "GlobalAdjoint uses the finite-difference saltation matrix at the event at " *
                    "t = $tf: ReverseDiff and ForwardDiff cannot differentiate the callback. " *
                    "Set `autojacvec = false` in the `sensealg` to use finite differences " *
                    "without this warning.", maxlog = 1
            )
            S, c = _event_saltation(
                callback, labels, tf, p, parent, P, sol, prob, dif, alg_idx;
                u_event = u_f, p_after, u_after
            )
            St, Sx = _saltation_maps(S)
        else
            St, Sx, c = maps
        end
        eb = copy(zero_e)
        eb[dif] .= c
        next === nothing && (eb .+= after_states[k] .- sol(tf; continuity = :right))
        eb[alg_idx] .= 0
        ea = first_at_time ? u_f .- sol(tf; continuity = :left) : copy(zero_e)
        (; a = tf, b = tf, ea, eb, St, Sx, paired = true, saves = saves[k])
    end
    return (;
        prob, sol, knots, jumps, windows = _merge_windows(windows; dif), calls = windows, dif,
        alg_idx, params,
    )
end

# The calls at the time of the saved entry `index` that act before it was saved, in order.
function _calls_before_entry(data, index)
    t = data.sol.t[index]
    return [w for w in data.calls if w.b == t && w.saves < index]
end

# The parameters of the forward solve on each interval between its event times (an affect
# can change them). The value before each call is in its snapshot; the value after the last
# event is the one of the integrator at the end. A change that a DAE reinitialization makes
# after an affect is also included, because the next snapshot comes after it.
function _parameter_record(calls, snapshots, p_end)
    times = unique(first.(calls))
    first_call = [findfirst(c -> first(c) == t, calls) for t in times]
    values = [[first(snapshots[k]) for k in first_call]; [p_end]]
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
function _unrecorded_data(prob, sol)
    dif, alg_idx = _semiexplicit_components(prob, "GlobalAdjoint")
    knots, jumps = _knot_jumps(sol)
    params = (; times = eltype(sol.t)[], values = [sol.prob.p])
    return (; prob, sol, knots, jumps, windows = Any[], calls = Any[], dif, alg_idx, params)
end
