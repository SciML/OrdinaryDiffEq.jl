```@meta
CollapsedDocStrings = true
```
# OrdinaryDiffEqLinear

Methods for semi-linear differential equations.

## Installation

To be able to access the solvers in `OrdinaryDiffEqLinear`, you must first install them use the Julia package manager:

```julia
using Pkg
Pkg.add("OrdinaryDiffEqLinear")
```

This will only install the solvers listed at the bottom of this page.
If you want to explore other solvers for your problem,
you will need to install some of the other libraries listed in the navigation bar on the left.

## Example usage

```julia
using OrdinaryDiffEqLinear, SciMLOperators
function update_func(A, u, p, t)
    A[1, 1] = 0
    A[2, 1] = sin(u[1])
    A[1, 2] = -1
    A[2, 2] = 0
end
A0 = ones(2, 2)
A = DiffEqArrayOperator(A0, update_func = update_func)
u0 = ones(2)
tspan = (0.0, 30.0)
prob = ODEProblem(A, u0, tspan)
sol = solve(prob, LieRK4(), dt = 1 / 4)
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
