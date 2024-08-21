# OrdinaryDiffEqLinear

Methods for semi-linear differential equations.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqLinear", "LieRK4")
```

## Full list of solvers

```@docs
MagnusMidpoint
MagnusLeapfrog
LieEuler
MagnusGauss4
MagnusNC6
MagnusGL6
MagnusGL8
MagnusNC8
MagnusGL4
RKMK2
RKMK4
LieRK4
CG2
CG3
CG4a
MagnusAdapt4
CayleyEuler
LinearExponential
```
