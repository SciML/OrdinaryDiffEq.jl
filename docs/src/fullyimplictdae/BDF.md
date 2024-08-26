```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqBDF

Multistep BDF methods, good for large stiff systems.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqBDF", "QNDF")
```

## Full list of solvers

### DAE

```@docs
DImplicitEuler
DABDF2
DFBDF
```
