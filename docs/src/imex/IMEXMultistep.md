```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqIMEXMultistep

Solvers if your system of ordinary differential equations can be split up into the sum of a stiff and non-stiff part.
These are "standard" low order IMEX methods common in many naive partial differential
equation applications.

## Installation

To be able to access the solvers in `OrdinaryDiffEqIMEXMultistep`, you must first install them use the Julia package manager:

```julia
using Pkg
Pkg.add("OrdinaryDiffEqIMEXMultistep")
```

This will only install the solvers listed at the bottom of this page.
If you want to explore other solvers for your problem,
you will need to install some of the other libraries listed in the navigation bar on the left.

## Example usage

```julia
using OrdinaryDiffEqIMEXMultistep
f1 = (u, p, t) -> 2u
f2 = (u, p, t) -> 2u
u0 = 1.0
tspan = (0.0, 1.0)
prob = SplitODEProblem(f1, f2, u0, tspan)
sol = solve(prob, CNAB2(), dt = 1 / 10)
```

## Full list of solvers

```@docs
CNAB2
CNLF2
```
