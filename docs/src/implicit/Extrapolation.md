# [OrdinaryDiffEqExtrapolation](@id StiffExtrapolation)

Solvers based on within method parallelism.
These solvers perform well for medium sized systems of ordinary differential equations, of about 20 to 500 equations,
at low tolerances.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqExtrapolation", "ImplicitEulerBarycentricExtrapolation")
```

## Full list of solvers

```@docs
ImplicitEulerExtrapolation
ImplicitDeuflhardExtrapolation
ImplicitHairerWannerExtrapolation
ImplicitEulerBarycentricExtrapolation
```
