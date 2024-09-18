```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqPRK

Explicit solvers optimized for a certain number of parallel calls of the system of ordinary differential equations `f`.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqPRK", "KuttaPRK2p5", fixed_timesteps = true)
```

## Full list of solvers

```@docs
KuttaPRK2p5
```
