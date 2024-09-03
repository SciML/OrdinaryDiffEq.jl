```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqQPRK

Explicit solvers optimized for a certain number of parallel calls of the system of ordinary differential equations `f`.
Particularly good at low tolerances, when using quad-precision arithmetic, `Float128`.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqQPRK", "QPRK98")
```

## Full list of solvers

```@docs
QPRK98
```
