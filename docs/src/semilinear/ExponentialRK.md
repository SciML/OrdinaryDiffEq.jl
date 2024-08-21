# OrdinaryDiffEqExponentialRK

Methods for semi-linear differential equations.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqExponentialRK", "Exprb43")
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
