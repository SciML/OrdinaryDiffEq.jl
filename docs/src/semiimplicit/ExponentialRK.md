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

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqExponentialRK", "Exprb43", fixed_timesteps = true)
```

## Full list of solvers

```@docs; canonical=false
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
