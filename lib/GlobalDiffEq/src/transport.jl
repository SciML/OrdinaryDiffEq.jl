"""
    GlobalErrorTransport(alg; gtol=nothing, transport_alg=alg,
                         autodiff=ADTypes.AutoForwardDiff(), maxiters=6,
                         safety=0.8, companion_abstol, companion_reltol)

Wrap an ODE algorithm with global error estimation and control based on the
linearized error-transport (first variational) equation. After a dense forward
solve producing the interpolant `P(t)`, the global error `ε(t) ≈ y(t) - P(t)`
is estimated by integrating the linear companion ODE

```math
ε'(t) = J(P(t), p, t) \\, ε(t) + d(t), \\qquad ε(t_0) = 0,
```

where `d(t) = f(P(t), p, t) - P'(t)` is the defect of the dense output and `J`
is the Jacobian of `f`, applied matrix-free through Jacobian-vector products
computed with `autodiff` (any `ADTypes` backend supported by
DifferentiationInterface). This is the error-transport approach of Shampine
(1986) and Berzins (1988), in the continuous defect-driven form summarized by
Lang and Verwer (2007).

Set the requested absolute endpoint error with `gtol`:

```julia
solve(prob, GlobalErrorTransport(Tsit5(); gtol = 1.0e-6))
```

The solver then tightens local tolerances until the estimated endpoint global
error 2-norm is at most `gtol`, like [`GlobalAdjoint`](@ref). Omit `gtol` when
using the algorithm only with [`global_error_estimate`](@ref). `transport_alg`
selects the solver for the companion equation, and `companion_abstol` /
`companion_reltol` control its tolerances.

This implementation supports forward-time, standard-mass-matrix ODEs with real
vector states and no callbacks. The wrapped solver must provide a dense
first-derivative interpolation.

## References

  - L. F. Shampine, Global error estimation with one-step methods, Computers &
    Mathematics with Applications 12A (1986).
  - M. Berzins, Global error estimation in the method of lines for parabolic
    equations, SIAM Journal on Scientific and Statistical Computing 9 (1988).
  - J. Lang and J. Verwer, On global error estimation and control for initial
    value problems, SIAM Journal on Scientific Computing 29 (2007).
"""
struct GlobalErrorTransport{A, TA, AD, G, V} <: GlobalDiffEqAlgorithm
    alg::A
    transport_alg::TA
    autodiff::AD
    gtol::G
    options::V
end

function GlobalErrorTransport(
        alg;
        transport_alg = alg,
        autodiff = ADTypes.AutoForwardDiff(),
        gtol = nothing,
        maxiters = 6,
        safety = 0.8,
        companion_abstol = gtol === nothing ? 1.0e-10 : min(gtol / 100, 1.0e-10),
        companion_reltol = gtol === nothing ? 1.0e-8 : min(gtol / 100, 1.0e-8)
    )
    options = _companion_options(
        gtol, maxiters, safety, companion_abstol, companion_reltol
    )
    return GlobalErrorTransport(alg, transport_alg, autodiff, gtol, options)
end

SciMLBase.allows_arbitrary_number_types(::GlobalErrorTransport) = false
SciMLBase.allowscomplex(::GlobalErrorTransport) = false
SciMLBase.isautodifferentiable(::GlobalErrorTransport) = false

function _transport_rhs(sol, autodiff)
    prob = sol.prob
    foop = _oop_rhs(prob)
    return let sol = sol, foop = foop, p = prob.p, backend = autodiff
        function (ε, _, t)
            u = sol(t, continuity = :right)
            du = sol(t, Val{1}, continuity = :right)
            defect = foop(u, p, t) - du
            jv = only(
                DifferentiationInterface.pushforward(
                    x -> foop(x, p, t), backend, u, (ε,)
                )
            )
            return jv + defect
        end
    end
end

"""
    global_error_estimate(prob, alg::GlobalErrorTransport; abstol=1e-6, reltol=1e-3, kwargs...)

Solve an ODE problem and estimate the 2-norm of its global error at the final
time by integrating the linearized error-transport equation driven by the
dense-output defect. `abstol` and `reltol` control the forward solve;
`companion_abstol` and `companion_reltol` control the companion error solve.
"""
function global_error_estimate(
        prob::SciMLBase.AbstractODEProblem, alg::GlobalErrorTransport, args...;
        abstol = 1.0e-6,
        reltol = 1.0e-3,
        companion_abstol = alg.options.companion_abstol,
        companion_reltol = alg.options.companion_reltol,
        kwargs...
    )
    return _companion_error_estimate(
        sol -> _transport_rhs(sol, alg.autodiff), "GlobalErrorTransport",
        prob, alg.alg, alg.transport_alg, args...;
        abstol, reltol, companion_abstol, companion_reltol, kwargs...
    )
end

function SciMLBase.__solve(
        prob::SciMLBase.AbstractODEProblem, alg::GlobalErrorTransport, args...;
        abstol = something(alg.gtol, 1.0e-6),
        reltol = something(alg.gtol, 1.0e-3),
        kwargs...
    )
    alg.gtol === nothing && throw(
        ArgumentError("GlobalErrorTransport requires a positive `gtol` constructor keyword")
    )
    _validate_tolerances(abstol, reltol, "local")
    haskey(kwargs, :callback) &&
        throw(ArgumentError("GlobalErrorTransport does not currently support callbacks"))
    estimator = (local_abstol, local_reltol) -> global_error_estimate(
        prob, alg, args...;
        abstol = local_abstol, reltol = local_reltol, kwargs...
    )
    return _refine_to_gtol(
        estimator, prob, alg.alg, alg.gtol, alg.options, args...;
        abstol, reltol, kwargs...
    )
end
