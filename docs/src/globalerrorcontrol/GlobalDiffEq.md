# GlobalDiffEq: Global Error Estimation and Control

Standard adaptive ODE solvers control the *local* error of each step. The
local tolerances only indirectly control the *global* (accumulated) error of
the solution, which can grow arbitrarily large over long integrations or on
unstable problems even when every step satisfies its local tolerance. The
GlobalDiffEq sublibrary provides solvers and solver wrappers that estimate
the global error, and in several cases control it to a requested global
tolerance.

To use these methods:

```julia
using GlobalDiffEq
```

[`GlobalRichardson`](@ref) wraps any fixed-step method in global Richardson
extrapolation over whole solves, interpreting `abstol` and `reltol` as global
tolerances. It is the most robust and most expensive option.

To *control* the endpoint global error to a tolerance `gtol`, wrap any
adaptive solver in [`GlobalAdjoint`](@ref) (adjoint-based, for endpoint
functionals; requires SciMLSensitivity and QuadGK to be loaded). It solves the
problem, estimates the endpoint global error, and tightens the local
tolerances until the requested global tolerance is met.

## Global error controlling wrappers

```@docs
GlobalRichardson
GlobalAdjoint
adjoint_error_estimate
```
