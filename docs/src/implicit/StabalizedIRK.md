# OrdinaryDiffEqStabalizedIRK

Explicit stabilized methods utilize an upper bound on the spectral radius of the Jacobian.
Users can supply an upper bound by specifying the keyword argument `eigen_est`, for example

```julia
`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`
```

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqStabalizedIRK", "IRKC")
```

## Full list of solvers

```@docs
IRKC
```

