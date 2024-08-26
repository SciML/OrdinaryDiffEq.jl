```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqHighOrderRK

Solvers for non-stiff problems at low tolerance.
However, the solvers in [`OrdinaryDiffEqVerner`](@ref OrdinaryDiffEqVerner) generally perform better at low tolerances.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqHighOrderRK", "DP8")
```

## Full list of solvers

```@docs
TanYam7
TsitPap8
DP8
PFRK87
```
