# OrdinaryDiffEqLinear

Methods for semi-linear differential equations.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqLinear", "LieRK4")
```

## Full list of solvers

### Time and State-Independent Solvers

```@docs
LinearExponential
```

### Time-Dependent and State-Independent Solvers

```@docs
MagnusMidpoint
MagnusLeapfrog
MagnusGauss4
MagnusNC6
MagnusGL6
MagnusGL8
MagnusNC8
MagnusGL4
```

### State-Dependent Solvers

```@docs
LieEuler
RKMK2
RKMK4
LieRK4
CG2
CG4a
MagnusAdapt4
CayleyEuler
```

### Time and State-Dependent Operators

```@docs
CG3
```
