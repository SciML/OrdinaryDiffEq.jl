# OrdinaryDiffEqFeagin

Preferred solvers for non-stiff problems at very low tolerance, `<1e-30`.
Best combined with preciser than `Float64` number types for the state, such as the `BigFloat` number type.
Note that the Feagin methods have a less robust error estimator than the Verner methods, and thus even for
very low tolerance problems the Verner methods (`Vern9`) may still be more efficient. In addition, at extremely
low tolerances the explicit extrapolation methods allow for arbitrarily high variable order stepping which will
also outperform the Feagin methods. As such, the Feagin methods may be useful in the Float128 precision range
but should be tested against other algorithms.

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
