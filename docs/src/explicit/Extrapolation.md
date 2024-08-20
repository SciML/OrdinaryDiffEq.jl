# OrdinaryDiffEqExtrapolation

Solvers based on within method parallelism.
The explicit solvers are outclassed by other explicit methods.
However, some [stiff extrapolation](@ref StiffExtrapolation) methods perform very well.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqExtrapolation", "ExtrapolationMidpointDeuflhard")
```

## Full list of solvers

```@docs
AitkenNeville
ExtrapolationMidpointDeuflhard
ExtrapolationMidpointHairerWanner
```