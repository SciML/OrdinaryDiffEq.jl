```@meta
CollapsedDocStrings = true
```
# OrdinaryDiffEqExtrapolation

Solvers based on within method parallelism, allowing multithreading of the solution across
different values of `f`.
The explicit extrapolation solvers are generally outclassed by other explicit methods.
However, some [stiff extrapolation](@ref StiffExtrapolation) methods perform very well if
the problem is sufficiently stiff.

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
