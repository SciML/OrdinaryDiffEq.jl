# OrdinaryDiffEqBDF

Multistep methods, good for large stiff systems.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqBDF", "QNDF")
```

## Full list of solvers

```@docs
ABDF2
QNDF
QNDF1
QNDF2
QBDF
QBDF1
QBDF2
MEBDF2
FBDF
```

### IMEX Multistep

```@docs
SBDF
SBDF2
SBDF3
SBDF4
```

### IMEX SDIRK

```@docs
IMEXEuler
IMEXEulerARK
```

### DAE

```@docs
DImplicitEuler
DABDF2
DFBDF
```
