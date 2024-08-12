# OrdinaryDiffEqRKN

Second order solvers.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqRKN", "Nystrom4")
```

## Full list of solvers

```@docs
IRKN3
IRKN4
Nystrom4
Nystrom4VelocityIndependent
Nystrom5VelocityIndependent
FineRKN4
FineRKN5
DPRKN4
DPRKN5
DPRKN6
DPRKN6FM
DPRKN8
DPRKN12
ERKN4
ERKN5
ERKN7
RKN4
```
