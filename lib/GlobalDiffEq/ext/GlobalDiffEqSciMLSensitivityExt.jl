module GlobalDiffEqSciMLSensitivityExt

import GlobalDiffEq, SciMLBase, SciMLSensitivity
import Accessors: @set
import LinearAlgebra: dot

# The adjoint error identity of GlobalAdjoint. For the dual -Mᵀ ż = Jᵀ z along the forward
# solution, z(T) = v,
#
#     ⟨v, M e(T)⟩ = z(t₀)ᵀ M e(t₀) + ∫ zᵀ r dt - Σⱼ z(tⱼ)ᵀ M [P]ⱼ + (event terms),
#
# with r = f(P) - M P' the defect of the dense output P and [P]ⱼ its jump at the knot tⱼ.
# The integral is the sensitivity of ⟨v, u(T)⟩ to a scalar forcing ε that adds r to the
# dynamics (`u' = f(u) + ε r(t)`): SciMLSensitivity returns it as the parameter gradient.

# Skip from the current time to `t`: drop the tstops in between (`first_tstop` gives tdir·t,
# and the integrator checks after the callbacks that it passed none), then move the time. A
# tstop at `t` itself stays; the caller handles a stop there before it skips.
function _skip_to!(integrator, t)
    while SciMLBase.has_tstop(integrator) &&
            SciMLBase.first_tstop(integrator) < integrator.tdir * t
        SciMLBase.pop_tstop!(integrator)
    end
    SciMLBase.set_t!(integrator, t)
    SciMLBase.u_modified!(integrator, true)
    return nothing
end

# The terminal value of the dual. An algebraic part of v moves to the differential part
# through the constraint at T: with A the Jacobian (or the secant Jacobian) of f at T,
# ⟨v, e⟩ = ⟨v_d - A_adᵀ w, e_d⟩ - wᵀ f_a(P(T⁻)), with w = A_aa⁻ᵀ v_a.
function _terminal_value(v, A, fa, dif, alg_idx)
    isempty(alg_idx) && return copy(v), zero(eltype(v))
    w = A[alg_idx, alg_idx]' \ v[alg_idx]
    vterm = zero(v)
    vterm[dif] .= v[dif] .- A[alg_idx, dif]' * w
    return vterm, -dot(w, fa[alg_idx])
end

# One backward solve of the dual along `forced_sol` (a solution object whose problem has the
# forced right-hand side), with the events of the forward solve. At each knot before the
# terminal time a callback records z and skips one ulp, so that the step below reads the
# forward state and P on its own side (and an FSAL solver does not reuse the evaluation on
# the knot). At the end b of each event window it records the pair term, applies the bridge
# z_d ← Sᵀ z_d, and moves to just before a: the dual does not enter the window. A stop at the
# time it moves to (a knot one ulp below another one) is handled in the same call. For a DAE
# the callback's initialization makes z_a consistent at the new time. Returns the integral,
# the knot and pair terms, and z(t₀).
function _bridged_adjoint(
        forced_sol, sensealg, adjoint_alg, vterm, terminal_time;
        M, dif, alg_idx, knots, jumps, windows, abstol, reltol, reinit_abstol
    )
    n = length(vterm)
    isdae = !isempty(alg_idx)
    terminal! = (out, u, p, t, i) -> (copyto!(out, vterm); nothing)
    aprob = SciMLSensitivity.ODEAdjointProblem(forced_sol, sensealg, adjoint_alg, [terminal_time], terminal!)
    inwindow(t) = any(w -> w.a <= t <= w.b, windows)
    outer_knots = filter(t -> t < terminal_time && !inwindow(t), knots)
    active = filter(w -> w.b < terminal_time, windows)
    stops = sort(unique([outer_knots; [w.b for w in active]]))
    z_at_knots = Dict{eltype(knots), typeof(vterm)}()
    pair_terms = eltype(vterm)[]
    # The knot record and the bridges at the stop t; returns the time to move to.
    cross! = function (integrator, t)
        target = nothing
        if insorted(t, outer_knots)
            z_at_knots[t] = integrator.u[1:n]
            target = prevfloat(t)
        end
        for w in active
            w.b == t || continue
            zb = integrator.u[dif]
            za = w.St(zb)
            push!(pair_terms, dot(zb, w.eb[dif]) - dot(za, w.ea[dif]))
            integrator.u[dif] .= za
            target = prevfloat(w.a)
        end
        return target
    end
    affect! = function (integrator)
        target = cross!(integrator, integrator.t)
        while target !== nothing && insorted(target, stops)
            target = cross!(integrator, target)
        end
        target === nothing || _skip_to!(integrator, target)
        return nothing
    end
    bridge = SciMLBase.DiscreteCallback(
        (u, t, integrator) -> insorted(t, stops), affect!;
        save_positions = (false, false),
        initializealg = isdae ? GlobalDiffEq.BrownFullBasicInit(abstol = reinit_abstol) : SciMLBase.NoInit()
    )
    asol = SciMLBase.solve(
        aprob, adjoint_alg; callback = bridge, tstops = stops,
        save_everystep = false, save_start = false, save_end = true, abstol, reltol
    )
    SciMLBase.successful_retcode(asol) || throw(
        ErrorException("the adjoint solve failed with retcode $(asol.retcode) at t = $(asol.t[end])")
    )
    integral = asol.u[end][end]
    knot_terms = [-dot(z_at_knots[t][dif], (M * J)[dif]) for (t, J) in zip(knots, jumps) if t in keys(z_at_knots)]
    return (; integral, knot_terms, pair_terms, z0 = asol.u[end][1:n])
end

# GlobalAdjoint: the dual of the error equation linearized at P, along the forward solution.
# Each forward event is a window of zero width with the saltation matrix of the event, so it
# has no pair term (the linearized limit). `data` is the recorded forward solve, and
# `terminal_index` the saved entry whose error is projected. An event time is saved twice;
# the second entry is the state after the event, so its error is e⁺ and the dual crosses the
# event first.
function GlobalDiffEq._adjoint_defect_projection(
        data, sensealg, adjoint_alg, direction;
        abstol, reltol, terminal_index = lastindex(data.sol.t)
    )
    (; prob, sol, knots, jumps, windows, dif, alg_idx, params) = data
    terminal_time = sol.t[terminal_index]
    sensealg === nothing ||
        sensealg isa SciMLSensitivity.InterpolatingAdjoint ||
        return _plain_adjoint_projection(data, sensealg, adjoint_alg, direction; abstol, reltol, terminal_time)
    # The default: a ReverseDiff VJP, and a numerical VJP (finite differences, which need no
    # AD through the model) if ReverseDiff cannot run the model, as `adjoint_sensitivities`
    # chooses for an out-of-place problem with callbacks.
    automatic = sensealg === nothing
    sensealg = something(sensealg, SciMLSensitivity.InterpolatingAdjoint(autojacvec = SciMLSensitivity.ReverseDiffVJP()))
    # The dual must follow the dense output P itself; a checkpointed adjoint solves the
    # forward problem again between its checkpoints, which gives another curve.
    sensealg.checkpointing && throw(
        ArgumentError(
            "GlobalAdjoint does not support InterpolatingAdjoint(checkpointing = true): the " *
                "adjoint must follow the stored forward solution"
        )
    )
    f = SciMLBase.unwrapped_f(prob.f)
    iip = SciMLBase.isinplace(prob)
    M = prob.f.mass_matrix
    T = sol.t[end]
    # SciMLSensitivity reads the forward solution with the right limit; r and the parameters
    # (which an affect can change) match it. At T there is no right limit (the derivative of
    # a Rosenbrock interpolant is Inf there).
    side(t) = t < T ? :right : :left
    p_at(t, s) = GlobalDiffEq._parameters_at(params, t, s)
    rhs(x, t, s = side(t)) = GlobalDiffEq._rhs(f, iip, x, p_at(t, s), t)
    r(t) = rhs(sol(t; continuity = side(t)), t) .- M * sol(t, Val{1}; continuity = side(t))
    F = (u, εp, t) -> rhs(u, t) .+ εp[1] .* r(t)
    forced = SciMLBase.ODEProblem(
        SciMLBase.ODEFunction{false}(F; mass_matrix = M), sol.prob.u0, sol.prob.tspan,
        [zero(eltype(sol.prob.u0))]
    )
    # The calls at the terminal time that act before the saved entry. With none, the entry is
    # the state before the events there (the left side). Otherwise the dual first crosses
    # those calls, in reverse order, with their source terms (as the bridge does), and the
    # entry is the state after them: after all calls at that time it is P(t⁺) up to the
    # saved state, else the state between two calls.
    before_entry = GlobalDiffEq._calls_before_entry(data, terminal_index)
    after_event = !isempty(before_entry)
    between = after_event && length(before_entry) < count(w -> w.b == terminal_time, data.calls)
    terminal_side = after_event ? :right : :left
    PT = between ? sol.u[terminal_index] : sol(terminal_time; continuity = terminal_side)
    A = isempty(alg_idx) ? nothing :
        GlobalDiffEq._rhs_jacobian(f, iip, PT, p_at(terminal_time, terminal_side), terminal_time)
    vterm, endpoint_term = _terminal_value(direction, A, rhs(PT, terminal_time, terminal_side), dif, alg_idx)
    # The identity gives the error against P(t∓); the saved state can differ from it (a dense
    # output that does not end at the step's state, such as Rosenbrock32's).
    shift = dot(direction, PT .- sol.u[terminal_index])
    for w in reverse(before_entry)
        za = w.St(vterm[dif])
        endpoint_term += dot(vterm[dif], w.eb[dif]) - dot(za, w.ea[dif])
        vterm[dif] .= za
    end
    # The dual starts one ulp before an interior terminal time, so that its first step reads
    # the forward solution and the parameters before the knot (as at the other knots).
    start = terminal_time < T ? prevfloat(terminal_time) : terminal_time
    solve_terms(sensealg, tol) = _bridged_adjoint(
        (@set sol.prob = forced), sensealg, adjoint_alg, vterm, start;
        M, dif, alg_idx, knots, jumps, windows, abstol = tol, reltol, reinit_abstol = abstol
    )
    # z(t₀)ᵀ M e(t₀): an initialization of a callback can change the state at t₀ without a
    # save, so that the dense output starts at another state.
    e0 = data.e0
    total(terms) = terms.integral + sum(terms.knot_terms; init = zero(terms.integral)) +
        sum(terms.pair_terms; init = zero(terms.integral)) + endpoint_term + shift +
        dot(terms.z0[dif], (M * e0)[dif])
    # The integral (last component of the adjoint state) is a quadrature of the size of the
    # error, which is not known before the solve. If the result is small against its abstol,
    # solve again with an abstol of 1e-6 times the result for the integral.
    value, sensealg = try
        total(solve_terms(sensealg, abstol)), sensealg
    catch err
        (automatic && !(err isa InterruptException)) || rethrow()
        @warn(
            "GlobalAdjoint: the adjoint solve with a ReverseDiff VJP failed; it uses a " *
                "numerical VJP (pass a `sensealg` to choose)", maxlog = 1
        )
        numerical = SciMLSensitivity.InterpolatingAdjoint(autojacvec = false, autodiff = false)
        total(solve_terms(numerical, abstol)), numerical
    end
    if abs(value) < 1.0e4 * abstol
        n = length(direction)
        value = total(solve_terms(sensealg, [fill(abstol, n); max(1.0e-6 * abs(value), floatmin(abs(value)))]))
    end
    return value
end

# Any other adjoint method: SciMLSensitivity's `adjoint_sensitivities`, for problems without
# events and mass matrices. It has no knot terms, so it is accurate only for a dense output
# that is continuous at the step ends.
function _plain_adjoint_projection(data, sensealg, adjoint_alg, direction; abstol, reltol, terminal_time)
    (; prob, sol, windows, alg_idx, params) = data
    (isempty(windows) && isempty(alg_idx)) || throw(
        ArgumentError(
            "GlobalAdjoint supports events and DAEs only with an InterpolatingAdjoint sensealg " *
                "(the default)"
        )
    )
    f_orig = SciMLBase.unwrapped_f(prob.f)
    realp = only(params.values)
    inplace = SciMLBase.isinplace(prob)
    forced_f = if inplace
        let f = f_orig, realp = realp, sol = sol
            function (du, u, forcing, t)
                f(du, u, realp, t)
                residual = similar(du)
                f(residual, sol(t), realp, t)
                deriv = sol(t, Val{1})
                @. du += forcing[1] * (residual - deriv)
                return nothing
            end
        end
    else
        let f = f_orig, realp = realp, sol = sol
            function (u, forcing, t)
                residual = f(sol(t), realp, t) .- sol(t, Val{1})
                return f(u, realp, t) .+ forcing[1] .* residual
            end
        end
    end
    forced_prob = SciMLBase.remake(
        sol.prob;
        f = SciMLBase.ODEFunction{inplace}(forced_f),
        p = [zero(eltype(sol.prob.u0))]
    )
    forced_sol = @set sol.prob = forced_prob
    terminal_gradient! = let direction = direction
        function (out, u, p, t, i)
            copyto!(out, direction)
            return nothing
        end
    end
    _, defect_gradient = SciMLSensitivity.adjoint_sensitivities(
        forced_sol, adjoint_alg;
        sensealg, t = [terminal_time], dgdu_discrete = terminal_gradient!, abstol, reltol
    )
    return only(defect_gradient)
end

end
