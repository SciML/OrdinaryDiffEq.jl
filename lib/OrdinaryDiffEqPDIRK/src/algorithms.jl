"""
    PDIRK44(; autodiff = AutoForwardDiff(), concrete_jac = nothing, linsolve = nothing,
        nlsolve = NLNewton(), extrapolant = :constant, threading = true)

Fourth-order, parallel diagonally implicit Runge-Kutta method. It evaluates the
independent diagonal stages concurrently when `threading = true`, making it useful
when stage solves are expensive and the problem can use CPU parallelism.

`PDIRK44` is non-adaptive. Supply `dt` when solving an `ODEProblem`.

# Fields

  - `linsolve`: optional linear solver configuration used by the implicit stages.
  - `nlsolve`: nonlinear solver algorithm used for the stage equations.
  - `extrapolant`: initial-stage extrapolation strategy.
  - `threading`: whether eligible internal stage work uses threaded execution.
  - `autodiff`: automatic-differentiation backend used to form Jacobians.
  - `concrete_jac`: controls whether a concrete Jacobian is cached when supported.

# Keywords

  - `autodiff = AutoForwardDiff()`: ADTypes backend for Jacobian construction.
  - `concrete_jac = nothing`: retain the default Jacobian-concretization policy. Set a
    boolean value to request or disable a concrete Jacobian where the problem supports it.
  - `linsolve = nothing`: linear solver algorithm or configuration. `nothing` selects
    the OrdinaryDiffEq default.
  - `nlsolve = NLNewton()`: nonlinear solver for each implicit stage.
  - `extrapolant = :constant`: extrapolation strategy used to initialize stage solves.
  - `threading = true`: enable threaded independent stage work. Set to `false` for
    serial execution.

# Examples

```julia
using OrdinaryDiffEqPDIRK
using SciMLBase: ODEProblem, solve

prob = ODEProblem((u, p, t) -> -10u, 1.0, (0.0, 1.0))
sol = solve(prob, PDIRK44(), dt = 0.1)
```

# References

A. Iserles and S. P. Norsett, "On the theory of parallel Runge-Kutta methods,"
*IMA Journal of Numerical Analysis* 10 (1990), pp. 463-488.
"""
struct PDIRK44{AD, F, F2, TO, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    threading::TO
    autodiff::AD
    concrete_jac::CJ
end
function PDIRK44(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant, threading = true
    )
    autodiff = _fixup_ad(autodiff)

    return PDIRK44(
        linsolve, nlsolve,
        extrapolant, threading, autodiff,
        _unwrap_val(concrete_jac)
    )
end
