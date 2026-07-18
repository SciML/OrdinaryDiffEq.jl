"""
    GlobalDefectCorrection(alg; gtol=nothing, correction_alg=alg, maxiters=6,
                           safety=0.8, companion_abstol, companion_reltol)

Wrap an ODE algorithm with global error estimation and control based on
solving for the correction (defect correction). After a dense forward solve
producing the interpolant `P(t)`, the global error `ε(t) ≈ y(t) - P(t)` is
estimated by integrating the nonlinear companion ODE

```math
ε'(t) = f(P(t) + ε(t), p, t) - P'(t), \\qquad ε(t_0) = 0,
```

which requires no Jacobian, only extra evaluations of `f` on perturbed
arguments. This is the "solving for the correction" technique of Zadunaisky
(1976) and Dormand, Duckers and Prince (1984), in the continuous per-step
dense-output form of Dormand, Lockyer, McGorrigan and Prince (1989); its
validity under standard local-error step control was proven by Calvo, Higham,
Montijano and Rández (1996).

Set the requested absolute endpoint error with `gtol`:

```julia
solve(prob, GlobalDefectCorrection(Tsit5(); gtol = 1.0e-6))
```

The solver then tightens local tolerances until the estimated endpoint global
error 2-norm is at most `gtol`, like [`GlobalAdjoint`](@ref). Omit `gtol` when
using the algorithm only with [`global_error_estimate`](@ref).
`correction_alg` selects the solver for the companion equation, and
`companion_abstol` / `companion_reltol` control its tolerances.

This implementation supports forward-time, standard-mass-matrix ODEs with real
vector states and no callbacks. The wrapped solver must provide a dense
first-derivative interpolation.

## References

  - P. E. Zadunaisky, On the estimation of errors propagated in the numerical
    integration of ordinary differential equations, Numerische Mathematik 27
    (1976).
  - J. R. Dormand, R. R. Duckers and P. J. Prince, Global error estimation with
    Runge-Kutta methods, IMA Journal of Numerical Analysis 4 (1984).
  - J. R. Dormand, M. A. Lockyer, N. E. McGorrigan and P. J. Prince, Global
    error estimation with Runge-Kutta triples, Computers & Mathematics with
    Applications 18 (1989).
"""
struct GlobalDefectCorrection{A, CA, G, V} <: GlobalDiffEqAlgorithm
    alg::A
    correction_alg::CA
    gtol::G
    options::V
end

function GlobalDefectCorrection(
        alg;
        correction_alg = alg,
        gtol = nothing,
        maxiters = 6,
        safety = 0.8,
        companion_abstol = gtol === nothing ? 1.0e-10 : min(gtol / 100, 1.0e-10),
        companion_reltol = gtol === nothing ? 1.0e-8 : min(gtol / 100, 1.0e-8)
    )
    options = _companion_options(
        gtol, maxiters, safety, companion_abstol, companion_reltol
    )
    return GlobalDefectCorrection(alg, correction_alg, gtol, options)
end

SciMLBase.allows_arbitrary_number_types(::GlobalDefectCorrection) = false
SciMLBase.allowscomplex(::GlobalDefectCorrection) = false
SciMLBase.isautodifferentiable(::GlobalDefectCorrection) = false

function _correction_rhs(sol)
    prob = sol.prob
    foop = _oop_rhs(prob)
    return let sol = sol, foop = foop, p = prob.p
        function (ε, _, t)
            u = sol(t, continuity = :right)
            du = sol(t, Val{1}, continuity = :right)
            return foop(u + ε, p, t) - du
        end
    end
end

"""
    global_error_estimate(prob, alg::GlobalDefectCorrection; abstol=1e-6, reltol=1e-3, kwargs...)

Solve an ODE problem and estimate the 2-norm of its global error at the final
time by solving for the correction of the dense numerical solution. `abstol`
and `reltol` control the forward solve; `companion_abstol` and
`companion_reltol` control the companion error solve.
"""
function global_error_estimate(
        prob::SciMLBase.AbstractODEProblem, alg::GlobalDefectCorrection, args...;
        abstol = 1.0e-6,
        reltol = 1.0e-3,
        companion_abstol = alg.options.companion_abstol,
        companion_reltol = alg.options.companion_reltol,
        kwargs...
    )
    return _companion_error_estimate(
        _correction_rhs, "GlobalDefectCorrection",
        prob, alg.alg, alg.correction_alg, args...;
        abstol, reltol, companion_abstol, companion_reltol, kwargs...
    )
end

function SciMLBase.__solve(
        prob::SciMLBase.AbstractODEProblem, alg::GlobalDefectCorrection, args...;
        abstol = something(alg.gtol, 1.0e-6),
        reltol = something(alg.gtol, 1.0e-3),
        kwargs...
    )
    alg.gtol === nothing && throw(
        ArgumentError("GlobalDefectCorrection requires a positive `gtol` constructor keyword")
    )
    _validate_tolerances(abstol, reltol, "local")
    haskey(kwargs, :callback) &&
        throw(ArgumentError("GlobalDefectCorrection does not currently support callbacks"))
    estimator = (local_abstol, local_reltol) -> global_error_estimate(
        prob, alg, args...;
        abstol = local_abstol, reltol = local_reltol, kwargs...
    )
    return _refine_to_gtol(
        estimator, prob, alg.alg, alg.gtol, alg.options, args...;
        abstol, reltol, kwargs...
    )
end
