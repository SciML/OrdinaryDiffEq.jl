module GlobalDiffEqSciMLSensitivityExt

import DiffEqCallbacks, GlobalDiffEq, LinearAlgebra, SciMLBase, SciMLSensitivity
import Accessors: @set
using SciMLLogging: @SciMLMessage
import LinearAlgebra: dot

# The adjoint methods that integrate the error integral with a quadrature callback.
const GaussAdjoints = Union{SciMLSensitivity.GaussAdjoint, SciMLSensitivity.GaussKronrodAdjoint}

# The adjoint method of the algebraic projection in the DAE saltation products
# (`GlobalDiffEq._algebraic_solver`).
projection_sensealg() =
    SciMLSensitivity.SteadyStateAdjoint(autojacvec = SciMLSensitivity.ReverseDiffVJP())

# The adjoint error identity of GlobalAdjoint. For the dual -Mᵀ ż = Jᵀ z along the forward
# solution, with z(T) = v,
#
#     ⟨v, M e(T)⟩ = z(t₀)ᵀ M e(t₀) + ∫ zᵀ r dt - Σⱼ z(tⱼ)ᵀ M [P]ⱼ + (event terms),
#
# where r = f(P) - M P' is the defect of the dense output P and [P]ⱼ its jump at knot tⱼ.
# The integral is the sensitivity of ⟨v, u(T)⟩ to a scalar forcing ε in
# `u' = f(u) + ε r(t)`, so SciMLSensitivity returns it as the parameter gradient.

# Move the integrator to `t`. The tstops before `t` are dropped first (`first_tstop` returns
# tdir·t), because the integrator checks after the callbacks that it passed no tstop. A
# tstop at `t` stays; the caller handles it before the skip.
function _skip_to!(integrator, t)
    while SciMLBase.has_tstop(integrator) &&
            SciMLBase.first_tstop(integrator) < integrator.tdir * t
        SciMLBase.pop_tstop!(integrator)
    end
    SciMLBase.set_t!(integrator, t)
    SciMLBase.u_modified!(integrator, true)
    return nothing
end

# The terminal value of the dual. The constraint at T moves an algebraic part of v to the
# differential part: with A the Jacobian (or secant Jacobian) of f at T and w = A_aa⁻ᵀ v_a,
# ⟨v, e⟩ = ⟨v_d - A_adᵀ w, e_d⟩ - wᵀ f_a(P(T⁻)).
function _terminal_value(v, A, fa, dif, alg_idx)
    isempty(alg_idx) && return copy(v), zero(eltype(v))
    w = A[alg_idx, alg_idx]' \ v[alg_idx]
    vterm = zero(v)
    vterm[dif] .= v[dif] .- A[alg_idx, dif]' * w
    return vterm, -dot(w, fa[alg_idx])
end

# One backward solve of the dual along `forced_sol` (whose problem has the forced right-hand
# side), with the events of the forward solve. With a GaussAdjoint or GaussKronrodAdjoint
# `sensealg`, the state is z only, and a quadrature callback integrates zᵀ r on each adjoint
# step. The adjoint steps stop at each knot, so zᵀ r is smooth on each step, and the fixed
# rule is accurate relative to the integral. The quadrature callback comes before the
# bridge, so it reads each step before the bridge changes z. With an InterpolatingAdjoint,
# the integral is the last component of the adjoint state, and the error control of the
# solve applies to it (an absolute accuracy of about the abstol). At each knot before the
# terminal time, a
# callback records z and skips one ulp. The next step then reads the forward state and P on
# its own side, and an FSAL solver does not reuse the evaluation at the knot. At the end b
# of each event window, the callback records the pair term, applies z_d ← Sᵀ z_d, and skips
# to just before a, so the dual never enters the window. If the new time is also a stop, the
# same call handles it. For a DAE, the callback initialization makes z_a consistent.
function _bridged_adjoint(
        forced_sol, sensealg, adjoint_alg, vterm, terminal_time;
        M, dif, alg_idx, knots, jumps, windows, abstol, reltol, reinit_abstol, defect
    )
    n = length(vterm)
    isdae = !isempty(alg_idx)
    terminal! = (out, u, p, t, i) -> (copyto!(out, vterm); nothing)
    gauss = sensealg isa GaussAdjoints
    if gauss
        # The integrand is -λᵀ r with r = ∂F/∂ε exact (the sign of SciMLSensitivity's
        # adjoint convention). SciMLSensitivity's own integrand (`GaussIntegrand`) gets
        # λᵀ ∂F/∂ε as a VJP; with `autojacvec = false` that is a finite difference in ε,
        # which loses the small r against the large f(u). The method of `ODEAdjointProblem`
        # for a Gauss adjoint still needs a `GaussIntegrand`.
        integrand = SciMLSensitivity.GaussIntegrand(forced_sol, sensealg, forced_sol.t, nothing)
        quad_values = DiffEqCallbacks.IntegrandValuesSum(zeros(eltype(vterm), 1))
        quad_f = (out, u, t, integrator) -> (out[1] = -dot(u, defect(t)); nothing)
        quad = sensealg isa SciMLSensitivity.GaussAdjoint ?
            DiffEqCallbacks.IntegratingSumCallback(quad_f, quad_values, zeros(eltype(vterm), 1); integrand_inplace = true) :
            DiffEqCallbacks.IntegratingGKSumCallback(quad_f, quad_values, zeros(eltype(vterm), 1); integrand_inplace = true)
        aprob, _, _ = SciMLSensitivity.ODEAdjointProblem(
            forced_sol, sensealg, adjoint_alg, integrand, quad, [terminal_time], terminal!,
            nothing, nothing, nothing, nothing, Val(false); abstol, reltol
        )
        extra = (quad,)
    else
        aprob = SciMLSensitivity.ODEAdjointProblem(forced_sol, sensealg, adjoint_alg, [terminal_time], terminal!)
        extra = ()
    end
    # The windows are sorted and do not overlap, so a binary search finds the window of a
    # time. A linear search costs O(N²) for a callback that acts at each step.
    starts = [w.a for w in windows]
    inwindow(t) = (k = searchsortedlast(starts, t); k > 0 && t <= windows[k].b)
    outer_knots = filter(t -> t < terminal_time && !inwindow(t), knots)
    active = filter(w -> w.b < terminal_time, windows)
    ends = [w.b for w in active]
    stops = sort(unique([outer_knots; ends]))
    z_at_knots = Dict{eltype(knots), typeof(vterm)}()
    pair_terms = eltype(vterm)[]
    # Handle the stop `t` (knot record and bridges); return the time to skip to, or nothing.
    cross! = function (integrator, t)
        target = nothing
        if insorted(t, outer_knots)
            z_at_knots[t] = integrator.u[1:n]
            target = prevfloat(t)
        end
        for k in searchsorted(ends, t)
            w = active[k]
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
        aprob, adjoint_alg; callback = SciMLBase.CallbackSet(extra..., bridge), tstops = stops,
        save_everystep = false, save_start = false, save_end = true, abstol, reltol
    )
    SciMLBase.successful_retcode(asol) || throw(
        ErrorException("the adjoint solve failed with retcode $(asol.retcode) at t = $(asol.t[end])")
    )
    integral = gauss ? quad_values.integrand[1] : asol.u[end][end]
    knot_terms = [-dot(z_at_knots[t][dif], (M * J)[dif]) for (t, J) in zip(knots, jumps) if t in keys(z_at_knots)]
    return (; integral, knot_terms, pair_terms, z0 = asol.u[end][1:n])
end

# The GlobalAdjoint projection of the error at the saved entry `terminal_index` of the
# recorded forward solve `data`. The dual is linearized at P. Each forward event is a window
# of zero width with its saltation matrix, so it has no pair term. An event time is saved
# twice; for the second entry (after the event), the dual crosses the event first.
function GlobalDiffEq._adjoint_defect_projection(
        data, sensealg, adjoint_alg, direction;
        abstol, reltol, terminal_index = lastindex(data.sol.t)
    )
    (; prob, sol, knots, jumps, windows, dif, alg_idx, params) = data
    terminal_time = sol.t[terminal_index]
    sensealg === nothing ||
        sensealg isa SciMLSensitivity.InterpolatingAdjoint ||
        sensealg isa GaussAdjoints ||
        return _plain_adjoint_projection(data, sensealg, adjoint_alg, direction; abstol, reltol, terminal_time)
    # The default is a ReverseDiff VJP, with a numerical VJP if ReverseDiff cannot run the
    # model. `adjoint_sensitivities` makes the same choice for an out-of-place problem with
    # callbacks.
    automatic = sensealg === nothing
    sensealg = something(sensealg, SciMLSensitivity.GaussKronrodAdjoint(autojacvec = SciMLSensitivity.ReverseDiffVJP()))
    # The dual must follow the dense output P. A checkpointed adjoint solves the forward
    # problem again between checkpoints and so follows another curve.
    sensealg.checkpointing && throw(
        ArgumentError(
            "GlobalAdjoint does not support an adjoint with `checkpointing = true`: the " *
                "adjoint must follow the stored forward solution"
        )
    )
    f = SciMLBase.unwrapped_f(prob.f)
    iip = SciMLBase.isinplace(prob)
    # `_semiexplicit_components` allows only the identity as a UniformScaling mass matrix,
    # and SciMLSensitivity recognizes the identity only as the default `I`.
    M = prob.f.mass_matrix isa LinearAlgebra.UniformScaling ? LinearAlgebra.I : prob.f.mass_matrix
    T = sol.t[end]
    # SciMLSensitivity reads the forward solution from the right, so r and the parameters
    # (an affect can change them) do too. At T there is no right limit: the derivative of a
    # Rosenbrock interpolant is Inf there.
    side(t) = t < T ? :right : :left
    p_at(t, s) = GlobalDiffEq._parameters_at(params, t, s)
    rhs(x, t, s = side(t)) = GlobalDiffEq._rhs(f, iip, x, p_at(t, s), t)
    r(t) = rhs(sol(t; continuity = side(t)), t) .- M * sol(t, Val{1}; continuity = side(t))
    F = (u, εp, t) -> rhs(u, t) .+ εp[1] .* r(t)
    forced = SciMLBase.ODEProblem(
        SciMLBase.ODEFunction{false}(F; mass_matrix = M), sol.prob.u0, sol.prob.tspan,
        [zero(eltype(sol.prob.u0))]
    )
    # Work around a SciMLSensitivity bug: the finite-difference Jacobian fails for an
    # out-of-place DAE (`autojacvec = false, autodiff = false` errors in `default_relstep`).
    # `autodiff` is not public, so every DAE with `autojacvec = false` uses an in-place view
    # of F; with `autodiff = true` the view gives the same result.
    F! = (du, u, εp, t) -> (du .= F(u, εp, t); nothing)
    forced_inplace = isempty(alg_idx) ? nothing : SciMLBase.ODEProblem(
            SciMLBase.ODEFunction{true}(F!; mass_matrix = M), forced.u0, forced.tspan, forced.p
        )
    # `before_entry` holds the event calls at the terminal time that act before the saved
    # entry. With none, the entry is the state on the left side. Otherwise the dual first
    # crosses these calls in reverse order, with their pair terms, as the bridge does. The
    # entry is then P(t⁺) up to the saved state, or the saved state between two calls.
    before_entry = GlobalDiffEq._calls_before_entry(data, terminal_index)
    after_event = !isempty(before_entry)
    between = after_event && length(before_entry) < count(w -> w.b == terminal_time, data.calls)
    terminal_side = after_event ? :right : :left
    PT = between ? sol.u[terminal_index] : sol(terminal_time; continuity = terminal_side)
    A = isempty(alg_idx) ? nothing :
        GlobalDiffEq._rhs_jacobian(f, iip, PT, p_at(terminal_time, terminal_side), terminal_time, data.verbose)
    vterm, endpoint_term = _terminal_value(direction, A, rhs(PT, terminal_time, terminal_side), dif, alg_idx)
    # The identity gives the error against P(t∓), and the saved state can differ from it
    # (for example, the dense output of Rosenbrock32 does not end at the step's state).
    shift = dot(direction, PT .- sol.u[terminal_index])
    for w in reverse(before_entry)
        bridge = GlobalDiffEq._entry_bridge(w, terminal_index)
        za = bridge.St(vterm[dif])
        endpoint_term += dot(vterm[dif], bridge.eb[dif]) - dot(za, w.ea[dif])
        vterm[dif] .= za
    end
    # The dual starts one ulp before an interior terminal time, so that its first step reads
    # the forward solution and the parameters before the knot (as at the other knots).
    start = terminal_time < T ? prevfloat(terminal_time) : terminal_time
    solve_with(problem, sensealg) = _bridged_adjoint(
        (@set sol.prob = problem), sensealg, adjoint_alg, vterm, start;
        M, dif, alg_idx, knots, jumps, windows, abstol, reltol, reinit_abstol = abstol,
        defect = r,
    )
    # One call per problem type, so that each adjoint solve has a concrete type.
    solve_terms(sensealg) = !isempty(alg_idx) && sensealg.autojacvec === false ?
        solve_with(forced_inplace, sensealg) : solve_with(forced, sensealg)
    # The term z(t₀)ᵀ M e(t₀) is nonzero if a callback initialization changes the state at
    # t₀ without a save, so that the dense output starts at another state.
    e0 = data.e0
    # A new binding: the closure would box `endpoint_term`, which is reassigned.
    fixed_terms = endpoint_term + shift
    total(terms) = terms.integral + sum(terms.knot_terms; init = zero(terms.integral)) +
        sum(terms.pair_terms; init = zero(terms.integral)) + fixed_terms +
        dot(terms.z0[dif], (M * e0)[dif])
    return try
        total(solve_terms(sensealg))
    catch err
        (automatic && !(err isa InterruptException)) || rethrow()
        @SciMLMessage(
            "GlobalAdjoint: the adjoint solve with a ReverseDiff VJP failed; it uses a " *
                "numerical VJP (pass a `sensealg` to choose)",
            data.verbose, :global_error_ad_fallback, maxlog = 1
        )
        numerical = SciMLSensitivity.GaussKronrodAdjoint(autojacvec = false, autodiff = false)
        # If the numerical VJP fails too, the first error is the likely cause (for example a
        # missing package for the DAE reinitialization). `rethrow()` throws it with its own
        # backtrace, and the second error goes to a warning.
        fallback = try
            total(solve_terms(numerical))
        catch second
            second isa InterruptException && rethrow()
            @SciMLMessage(
                "GlobalAdjoint: the adjoint solve with a numerical VJP failed too; the error " *
                    "of the ReverseDiff VJP follows",
                data.verbose, :global_error_ad_fallback, exception = (second, catch_backtrace())
            )
            nothing
        end
        fallback === nothing && rethrow()
        fallback
    end
end

GlobalDiffEq._supports_events(::SciMLSensitivity.InterpolatingAdjoint) = true
GlobalDiffEq._supports_events(::GaussAdjoints) = true

# Any other adjoint method uses SciMLSensitivity's `adjoint_sensitivities`, without events
# or mass matrices. It has no knot terms, so it is accurate only for a dense output that is
# continuous at the step ends. `_recorded_forward_solve` rejects events and DAEs for it; a
# call that changes nothing, such as a SavingCallback, is no event. The fixed-grid solution
# of StepGridRefinement has no events but can be a DAE, so the check below rejects it.
function _plain_adjoint_projection(data, sensealg, adjoint_alg, direction; abstol, reltol, terminal_time)
    (; prob, sol, alg_idx, params) = data
    isempty(alg_idx) || throw(GlobalDiffEq._unsupported_events_error())
    f_orig = SciMLBase.unwrapped_f(prob.f)
    # No call changed the parameters, so all values in the record are equal.
    realp = first(params.values)
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
