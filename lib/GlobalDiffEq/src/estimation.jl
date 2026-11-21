"""
    GlobalErrorEstimation(alg; equation=DefectCorrection, mode=InterpolatingMode,
                          companion_alg=alg, autodiff=ADTypes.AutoForwardDiff(),
                          gtol=nothing, maxiters=6, safety=0.8,
                          companion_abstol, companion_reltol)

Wrap an ODE algorithm with global error estimation and control based on a
defect-driven companion equation. After a dense forward solve producing the
interpolant `P(t)`, the global error `ε(t) ≈ y(t) - P(t)` is estimated by
integrating a companion ODE driven by the dense-output defect
`d(t) = f(P(t), p, t) - P'(t)`, with `ε(t₀) = 0`.

Two orthogonal choices control the estimator:

`equation` ([`GlobalErrorEquation`](@ref)) selects the companion ODE:

  - [`DefectCorrection`](@ref) (default): the nonlinear correction
    `ε' = f(P + ε, p, t) - P'(t)`, requiring no Jacobian.
  - [`ErrorTransport`](@ref): the linearized transport `ε' = J(P, p, t) ε + d(t)`,
    with `J = ∂f/∂u` applied matrix-free as a Jacobian-vector product through
    `autodiff` (any `ADTypes` backend supported by DifferentiationInterface).
    `ErrorTransport` is the first-order linearization of `DefectCorrection`.

`mode` ([`GlobalErrorMode`](@ref)) selects how the companion is integrated:

  - [`InterpolatingMode`](@ref) (default): run the forward solve with dense
    output and `save_everystep`, keep it, and integrate the companion against its
    interpolant in a second solve. Solver-agnostic but `O(n_steps)` memory.
  - [`SimultaneousMode`](@ref): drive the forward integrator one accepted step at
    a time and advance the companion across each step using that step's live dense
    output, storing no dense solution (`O(1)` extra memory) and running no second
    solve. Its estimate closely agrees with `InterpolatingMode`.

Set the requested absolute endpoint error with `gtol`:

```julia
solve(prob, GlobalErrorEstimation(Tsit5(); gtol = 1.0e-6))
```

The solver then tightens local tolerances until the estimated endpoint global
error 2-norm is at most `gtol`. Without `gtol`, the solve does one estimate at the
given `abstol` and `reltol`, with no control.

With `equation = DefectCorrection`, the returned solution has the error trajectory
(`SciMLBase.has_global_error` is `true`): `sol.global_error[i]` is the estimated global
error `y(sol.t[i]) - sol.u[i]` of each saved state, as a vector. The solution is dense
and saved at each step of the forward solve (`saveat` does not apply). If an event time
is saved twice, the first entry has the error before the events at that time and the
second the error after them. Under global error control, the final solution costs one
more forward solve in `InterpolatingMode` and one more estimate in `SimultaneousMode`.
If the solve is deterministic, the user's callback objects see only that final solve.
[`global_error_estimate`](@ref) gives only the endpoint error 2-norm, with no solution.

`companion_alg` selects the solver for the companion equation (by default the forward
solver `alg`), and `companion_abstol` / `companion_reltol` control its tolerances. The
companion runs at much tighter tolerances than the forward solve. So with a forward solver
of order 3 or lower (such as Rosenbrock23 or Rosenbrock32), pass a higher-order
`companion_alg`, for example `Tsit5()`, or `Rodas5P()` for a stiff problem or a DAE. The
solve warns about such a low-order companion solver. The companion state is the global
error itself. `companion_reltol` defaults to `1e-8`, because the error can have a
transient peak much larger than its final value (for example at a fast transition of a
stiff problem). `companion_abstol` must be below the size of the error. It defaults to
`1e-4` times the forward `abstol`, and at most `gtol / 100` under global error control.

Every combination **requires** the wrapped solver to provide a genuine dense
first-derivative interpolation (a solver with only the trivial fallback
interpolation is rejected), a forward time span, and a real vector state.

`DefectCorrection` (the default, in either mode) also supports mass-matrix DAEs and
callbacks (`ContinuousCallback`, `VectorContinuousCallback`, `DiscreteCallback`), for
example ModelingToolkit models with events. The companion applies the user's callbacks
to the true state `P + ε`, with its own copy of the parameters. So the estimate includes
the effect of the event times and of events that the forward solve missed or added. A
warning reports when the events of the forward solve and of the true solution differ. A
`DiscreteCallback` is checked only at the step ends of the forward solve and at the times
of the true continuous events.

After an event in a DAE, the companion reinitializes the true state with the algorithm of
the forward solve (for example `BrownFullBasicInit`). `OverrideInit`, the default for a
ModelingToolkit problem when a callback has no `initializealg`, is not supported there.
A callback can terminate the solve; the estimate is then the error at the time where the
forward solve stops. A warning reports when the true solution does not terminate at the
same event, or terminates where the forward solve does not. Callbacks that change the
step size are not supported in the companion. `ErrorTransport` supports only
standard-mass-matrix ODEs without callbacks.

Warnings follow the `verbose` keyword of `solve` and of `global_error_estimate` (a
[`DEVerbosity`](@ref)), with these toggles: `global_error_event_mismatch` (the events or
the terminations of the forward solve and of the true solution differ),
`global_error_companion_order` (a companion solver of order 3 or lower) and
`global_error_final_solve` (the final forward solve does not repeat the accepted one).

## References

  - P. E. Zadunaisky, On the estimation of errors propagated in the numerical
    integration of ordinary differential equations, Numerische Mathematik 27
    (1976).
  - J. R. Dormand, R. R. Duckers and P. J. Prince, Global error estimation with
    Runge-Kutta methods, IMA Journal of Numerical Analysis 4 (1984).
  - L. F. Shampine, Global error estimation with one-step methods, Computers &
    Mathematics with Applications 12A (1986).
  - J. Lang and J. Verwer, On global error estimation and control for initial
    value problems, SIAM Journal on Scientific Computing 29 (2007).
"""
struct GlobalErrorEstimation{A, CA, AD, G, V} <: GlobalDiffEqAlgorithm
    alg::A
    companion_alg::CA
    equation::GlobalErrorEquation
    mode::GlobalErrorMode
    autodiff::AD
    gtol::G
    options::V
end

function GlobalErrorEstimation(
        alg;
        equation = DefectCorrection,
        mode = InterpolatingMode,
        companion_alg = nothing,
        autodiff = ADTypes.AutoForwardDiff(),
        gtol = nothing,
        maxiters = 6,
        safety = 0.8,
        companion_abstol = nothing,
        companion_reltol = nothing
    )
    (equation === DefectCorrection || equation === ErrorTransport) || throw(
        ArgumentError("`equation` must be DefectCorrection or ErrorTransport")
    )
    (mode === InterpolatingMode || mode === SimultaneousMode) || throw(
        ArgumentError("`mode` must be InterpolatingMode or SimultaneousMode")
    )
    options = _companion_options(
        gtol, maxiters, safety, companion_abstol, companion_reltol
    )
    return GlobalErrorEstimation(
        alg, something(companion_alg, alg), equation, mode, autodiff, gtol, options
    )
end

# The companion runs at much tighter tolerances than the forward problem, and a low-order
# solver meets them poorly or with very many steps.
function _warn_if_low_order_companion(alg::GlobalErrorEstimation, kwargs)
    order = try
        SciMLBase.alg_order(alg.companion_alg)
    catch
        return nothing
    end
    order isa Integer && order <= 3 && @SciMLMessage(
        "GlobalErrorEstimation solves the companion with $(nameof(typeof(alg.companion_alg))) " *
            "(order $order) at tight tolerances, which can give a poor estimate. Pass a " *
            "higher-order `companion_alg`, for example `Tsit5()`, or `Rodas5P()` for a stiff " *
            "problem or a DAE.",
        _verbosity(kwargs), :global_error_companion_order
    )
    return nothing
end

SciMLBase.allows_arbitrary_number_types(::GlobalErrorEstimation) = false
SciMLBase.allowscomplex(::GlobalErrorEstimation) = false
SciMLBase.isautodifferentiable(::GlobalErrorEstimation) = false
# The outer solve wraps `f` based on these traits, and the inner solves keep that wrapper.
# Without them, an implicit solver finds no dual signature.
SciMLBase.forwarddiffs_model(alg::GlobalErrorEstimation) = SciMLBase.forwarddiffs_model(alg.alg)
SciMLBase.forwarddiff_chunksize(alg::GlobalErrorEstimation) = SciMLBase.forwarddiff_chunksize(alg.alg)

"""
    global_error_estimate(prob, alg::GlobalErrorEstimation; abstol=1e-6, reltol=1e-3, kwargs...)

Solve an ODE problem and estimate the 2-norm of its global error at the final
time using the companion equation selected by `alg.equation` and integrated in
the mode selected by `alg.mode`. `abstol` and `reltol` control the forward solve;
`companion_abstol` and `companion_reltol` control the companion error solve.
"""
function global_error_estimate(
        prob::SciMLBase.AbstractODEProblem, alg::GlobalErrorEstimation, args...;
        abstol = 1.0e-6,
        reltol = 1.0e-3,
        companion_abstol = alg.options.companion_abstol,
        companion_reltol = alg.options.companion_reltol,
        kwargs...
    )
    _warn_if_low_order_companion(alg, kwargs)
    return _estimate(prob, alg, args...; abstol, reltol, companion_abstol, companion_reltol, kwargs...)
end

SciMLBase.has_global_error(alg::GlobalErrorEstimation) = alg.equation === DefectCorrection

# With `trajectory = true` (DefectCorrection only), returns `(; estimate, sol)`, where `sol`
# has the estimated error of each saved state in `global_error`.
function _estimate(
        prob, alg::GlobalErrorEstimation, args...;
        abstol, reltol, companion_abstol = alg.options.companion_abstol,
        companion_reltol = alg.options.companion_reltol, trajectory = false, kwargs...
    )
    companion_abstol, companion_reltol = _companion_tolerances(
        companion_abstol, companion_reltol, abstol, alg.gtol
    )
    simultaneous = alg.mode === SimultaneousMode
    if alg.equation === DefectCorrection
        estimate = simultaneous ? _defect_correction_estimate_streaming : _defect_correction_estimate
        return estimate(
            "GlobalErrorEstimation", prob, alg.alg, alg.companion_alg, args...;
            abstol, reltol, companion_abstol, companion_reltol, trajectory, kwargs...
        )
    end
    make_rhs = _companion_rhs_builder(alg.equation, alg.autodiff)
    estimate = simultaneous ? _companion_error_estimate_streaming : _companion_error_estimate
    return estimate(
        make_rhs, "GlobalErrorEstimation",
        prob, alg.alg, alg.companion_alg, args...;
        abstol, reltol, companion_abstol, companion_reltol, kwargs...
    )
end

# The returned solution under global error control with DefectCorrection. The accepted
# estimate ran on a copy of the callbacks, so this solve runs on the user's callbacks. In
# InterpolatingMode, only the forward solve runs again and reuses the accepted error
# trajectory, because a companion solve costs several forward solves. In SimultaneousMode,
# `accepted` is only the norm, so the full estimate runs again.
function _final_solve(prob, alg, accepted, args...; abstol, reltol, kwargs...)
    full_estimate(kw) = _estimate(prob, alg, args...; abstol, reltol, trajectory = true, kw...).sol
    accepted isa Number && return full_estimate(kwargs)
    # Copy the callbacks before the forward solve changes them, for the fallback below.
    unused = _with_copied_callback(kwargs)
    # As in the estimate, the forward affects change their own copy of the parameters.
    sol = SciMLBase.solve(
        (@set prob.p = deepcopy(prob.p)), alg.alg, args...;
        abstol, reltol, merge((; kwargs...), _DENSE_SOLVE_KWARGS)...
    )
    sol.t == accepted.sol.t && sol.u == accepted.sol.u &&
        return @set sol.global_error = accepted.sol.global_error
    @SciMLMessage(
        "GlobalErrorEstimation: the final forward solve does not repeat the accepted one (a " *
            "callback with state outside of it, or a solve that is not deterministic). The " *
            "returned solution comes from a new estimate on a copy of the callbacks, so the " *
            "callback objects saw another solve.",
        _verbosity(kwargs), :global_error_final_solve
    )
    return full_estimate(unused)
end

function SciMLBase.__solve(
        prob::SciMLBase.AbstractODEProblem, alg::GlobalErrorEstimation, args...;
        abstol = something(alg.gtol, 1.0e-6),
        reltol = something(alg.gtol, 1.0e-3),
        kwargs...
    )
    alg.gtol === nothing && alg.equation !== DefectCorrection && throw(
        ArgumentError(
            "GlobalErrorEstimation with equation = ErrorTransport requires a positive " *
                "`gtol` constructor keyword"
        )
    )
    _validate_tolerances(abstol, reltol, "local")
    alg.equation === DefectCorrection || !_has_callbacks(kwargs) ||
        throw(_callback_support_error("GlobalErrorEstimation"))
    _warn_if_low_order_companion(alg, kwargs)
    # The solve pipeline has already merged the problem's callback into `kwargs`. Remove it
    # from the problem, or the inner solves merge it again and an affect runs twice.
    prob = @set prob.kwargs = _without_callback(prob.kwargs)
    alg.gtol === nothing &&
        return _estimate(prob, alg, args...; abstol, reltol, trajectory = true, kwargs...).sol
    # Each estimate runs on a copy of the callbacks, so that the user's callback objects
    # (for example the values of a SavingCallback) see only the final solve. Only
    # InterpolatingMode keeps the trajectory; SimultaneousMode keeps its O(1) memory.
    trajectory = alg.equation === DefectCorrection && alg.mode === InterpolatingMode
    estimator = (local_abstol, local_reltol) -> _estimate(
        prob, alg, args...;
        abstol = local_abstol, reltol = local_reltol, trajectory, _with_copied_callback(kwargs)...
    )
    final_solve = alg.equation === DefectCorrection ?
        (local_abstol, local_reltol, accepted) -> _final_solve(
            prob, alg, accepted, args...; abstol = local_abstol, reltol = local_reltol, kwargs...
        ) : nothing
    return _refine_to_gtol(
        estimator, prob, alg.alg, alg.gtol, alg.options, args...;
        final_solve, abstol, reltol, kwargs...
    )
end
