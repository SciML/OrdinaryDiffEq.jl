# OrdinaryDiffEqFeagin

Preferred solvers for non-stiff problems at very low tolerance, `<1e-30`.
Best combined with preciser than `Float64` number types for the state, such as the `BigFloat` number type. 

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqFeagin", "Feagin14")
```

## Full list of solvers

```@docs
Feagin10
Feagin12
Feagin14
```