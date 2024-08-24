```@meta
CollapsedDocStrings = true
```
# [OrdinaryDiffEqVerner](@id OrdinaryDiffEqVerner)

Preferred solvers for non-stiff problems at low tolerance.
`Vern6`, `Vern7`, or `Vern8` are good methods for tolerances between `~1e-8-1e-12`,
and using `Float64` numbers for the state of the differential equation.
For even lower tolerances,`Vern9` should be used, combined with the more precise `BigFloat` number type.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqVerner", "Vern6")
```

## Full list of solvers

```@docs
Vern6
Vern7
Vern8
Vern9
```

```@docs
AutoVern6
AutoVern7
AutoVern8
AutoVern9
```
