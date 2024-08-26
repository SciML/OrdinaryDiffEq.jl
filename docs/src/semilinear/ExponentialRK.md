```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqExponentialRK

Methods for semi-linear differential equations.

## Installation

To be able to access the solvers in `OrdinaryDiffEqLinear`, you must first install them use the Julia package manager:

```julia
using Pkg
Pkg.add("OrdinaryDiffEqExponentialRK")
```

This will only install the solvers listed at the bottom of this page.
If you want to explore other solvers for your problem,
you will need to install some of the other libraries listed in the navigation bar on the left.

## Example usage

```julia
using OrdinaryDiffEqExponentialRK, SciMLOperators
A = [2.0 -1.0; -1.0 2.0]
linnonlin_f1 = MatrixOperator(A)
linnonlin_f2 = (du, u, p, t) -> du .= 1.01 .* u
linnonlin_fun_iip = SplitFunction(linnonlin_f1, linnonlin_f2)
tspan = (0.0, 1.0)
u0 = [0.1, 0.1]
prob = SplitODEProblem(linnonlin_fun_iip, u0, tspan)
sol = solve(prob, ETDRK4(), dt = 1 / 4)
```

## Full list of solvers

```@docs
LawsonEuler
NorsettEuler
ETD2
ETDRK2
ETDRK3
ETDRK4
HochOst4
```

### Adaptive Exponential Rosenbrock Methods

```@docs
Exprb32
Exprb43
```

### Exponential Propagation Iterative Runge-Kutta Methods (EPIRK)

```@docs
Exp4
EPIRK4s3A
EPIRK4s3B
EPIRK5s3
EXPRB53s3
EPIRK5P1
EPIRK5P2
```
