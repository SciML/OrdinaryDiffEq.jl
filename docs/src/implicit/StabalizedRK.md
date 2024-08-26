```@meta
CollapsedDocStrings = true
```
# OrdinaryDiffEqStabalizedRK

Explicit stabilized methods utilize an upper bound on the spectral radius of the Jacobian.
Users can supply an upper bound by specifying the keyword argument `eigen_est`, for example

```julia
`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`
```

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqStabalizedRK", "ROCK4")
```

## Full list of solvers

```@docs
ROCK2 
ROCK4 
RKC
SERK2
ESERK4
ESERK5
```
