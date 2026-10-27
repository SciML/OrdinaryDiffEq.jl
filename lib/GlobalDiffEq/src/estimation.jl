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
(`SciMLBase.has_global_error` is `true`): `sol.global_error[i]` is the estimated
global error `y(sol.t[i]) - sol.u[i]` of each saved state, as a vector. The solution is
dense and saved at each step of the forward solve (`saveat` does not apply). Where it
saves an event time two times, the first entry has the error before the events at that
time, and the second entry the error after them. Under global error control, the final
solve that gives the trajectory costs one more companion solve.
[`global_error_estimate`](@ref) gives only the endpoint error 2-norm, with no
solution. `companion_alg` selects the solver for the
companion equation, and `companion_abstol` / `companion_reltol` control its
tolerances. The companion state is the global error itself. `companion_reltol`
defaults to `1e-8`: the error can pass through a transient peak much larger than its
final value (for example at a fast transition of a stiff problem), so a loose relative
tolerance can give a wrong estimate. `companion_abstol` must be below the size of the
error and defaults to `1e-4` times the forward `abstol` (and at most `gtol / 100` under
global error control).

Every combination **requires** the wrapped solver to provide a genuine dense
first-derivative interpolation (a solver with only the trivial fallback
interpolation is rejected), a forward time span, and a real vector state.

`DefectCorrection` (the default, in either mode) also supports mass-matrix DAEs and
callbacks (`ContinuousCallback`, `VectorContinuousCallback`,
`DiscreteCallback`), for example ModelingToolkit models with events. The companion
then applies the user's callbacks to the true state `P + ε` with its own copy of the
parameters, so the estimate includes the effect of the event times and of events that
the forward solve missed or added; a warning reports when the events of the forward
solve and of the true solution differ. A `DiscreteCallback` is checked only at the step ends of the forward solve.
After an event in a DAE, the companion reinitializes the true state with the same
algorithm as the forward solve (for example `BrownFullBasicInit`); `OverrideInit`, the
default for a ModelingToolkit problem when a callback has no `initializealg`, is not
supported there. Callbacks that terminate the solve or change its step size are not
supported in the companion. `ErrorTransport` supports only standard-mass-matrix ODEs without callbacks.

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
        companion_alg = alg,
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
    return GlobalErrorEstimation(alg, companion_alg, equation, mode, autodiff, gtol, options)
end

SciMLBase.allows_arbitrary_number_types(::GlobalErrorEstimation) = false
SciMLBase.allowscomplex(::GlobalErrorEstimation) = false
SciMLBase.isautodifferentiable(::GlobalErrorEstimation) = false
# The outer solve wraps `f` for these traits, and the inner solves of the wrapped
# algorithm keep that wrapper; without them an implicit solver finds no dual signature.
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
    companion_abstol, companion_reltol = _companion_tolerances(
        companion_abstol, companion_reltol, abstol, alg.gtol
    )
    if alg.equation === DefectCorrection
        estimate = alg.mode === SimultaneousMode ?
            _defect_correction_estimate_streaming : _defect_correction_estimate
        return estimate(
            "GlobalErrorEstimation", prob, alg.alg, alg.companion_alg, args...;
            abstol, reltol, companion_abstol, companion_reltol, kwargs...
        )
    end
    make_rhs = _companion_rhs_builder(alg.equation, alg.autodiff)
    estimate = alg.mode === SimultaneousMode ?
        _companion_error_estimate_streaming : _companion_error_estimate
    return estimate(
        make_rhs, "GlobalErrorEstimation",
        prob, alg.alg, alg.companion_alg, args...;
        abstol, reltol, companion_abstol, companion_reltol, kwargs...
    )
end

SciMLBase.has_global_error(alg::GlobalErrorEstimation) = alg.equation === DefectCorrection

# The forward solution with the estimated error of each saved state in `global_error`
# (DefectCorrection only).
function _solve_with_global_error(prob, alg::GlobalErrorEstimation, args...; abstol, reltol, kwargs...)
    companion_abstol, companion_reltol = _companion_tolerances(
        alg.options.companion_abstol, alg.options.companion_reltol, abstol, alg.gtol
    )
    estimate = alg.mode === SimultaneousMode ?
        _defect_correction_estimate_streaming : _defect_correction_estimate
    return estimate(
        "GlobalErrorEstimation", prob, alg.alg, alg.companion_alg, args...;
        abstol, reltol, companion_abstol, companion_reltol, trajectory = true, kwargs...
    )
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
    # The solve pipeline has already merged the problem's callback into `kwargs`; remove
    # it from the problem so that the inner solves do not merge it a second time (an
    # affect such as x ← 2x would then run twice).
    prob = @set prob.kwargs = _without_callback(prob.kwargs)
    # Without gtol: one estimate at the given tolerances, with no control.
    alg.gtol === nothing &&
        return _solve_with_global_error(prob, alg, args...; abstol, reltol, kwargs...)
    # Each estimate runs on a copy of the callbacks, so that the user's callback objects
    # (for example the values of a SavingCallback) see only the final solve.
    estimator = (local_abstol, local_reltol) -> global_error_estimate(
        prob, alg, args...;
        abstol = local_abstol, reltol = local_reltol, _with_copied_callback(kwargs)...
    )
    # The final solve gives the error trajectory of the accepted solution. This costs one
    # more companion solve than a plain final solve.
    final_solve = alg.equation === DefectCorrection ?
        (local_abstol, local_reltol) -> _solve_with_global_error(
            prob, alg, args...; abstol = local_abstol, reltol = local_reltol, kwargs...
        ) : nothing
    return _refine_to_gtol(
        estimator, prob, alg.alg, alg.gtol, alg.options, args...;
        final_solve, abstol, reltol, kwargs...
    )
end
