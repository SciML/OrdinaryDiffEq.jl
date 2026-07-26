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

## Global error controlling wrappers

```@docs
GlobalRichardson
```
