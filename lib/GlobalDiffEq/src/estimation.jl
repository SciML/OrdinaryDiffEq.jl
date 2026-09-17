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
error 2-norm is at most `gtol`. Omit `gtol` when using the algorithm only with
[`global_error_estimate`](@ref). `companion_alg` selects the solver for the
companion equation, and `companion_abstol` / `companion_reltol` control its
tolerances.

Every combination **requires** the wrapped solver to provide a genuine dense
first-derivative interpolation (a solver with only the trivial fallback
interpolation is rejected). This implementation supports forward-time,
standard-mass-matrix ODEs with real vector states and no callbacks.

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
        companion_abstol = gtol === nothing ? 1.0e-10 : min(gtol / 100, 1.0e-10),
        companion_reltol = gtol === nothing ? 1.0e-8 : min(gtol / 100, 1.0e-8)
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
    make_rhs = _companion_rhs_builder(alg.equation, alg.autodiff)
    estimate = alg.mode === SimultaneousMode ?
        _companion_error_estimate_streaming : _companion_error_estimate
    return estimate(
        make_rhs, "GlobalErrorEstimation",
        prob, alg.alg, alg.companion_alg, args...;
        abstol, reltol, companion_abstol, companion_reltol, kwargs...
    )
end

function SciMLBase.__solve(
        prob::SciMLBase.AbstractODEProblem, alg::GlobalErrorEstimation, args...;
        abstol = something(alg.gtol, 1.0e-6),
        reltol = something(alg.gtol, 1.0e-3),
        kwargs...
    )
    alg.gtol === nothing && throw(
        ArgumentError("GlobalErrorEstimation requires a positive `gtol` constructor keyword")
    )
    _validate_tolerances(abstol, reltol, "local")
    haskey(kwargs, :callback) &&
        throw(ArgumentError("GlobalErrorEstimation does not currently support callbacks"))
    estimator = (local_abstol, local_reltol) -> global_error_estimate(
        prob, alg, args...;
        abstol = local_abstol, reltol = local_reltol, kwargs...
    )
    return _refine_to_gtol(
        estimator, prob, alg.alg, alg.gtol, alg.options, args...;
        abstol, reltol, kwargs...
    )
end
