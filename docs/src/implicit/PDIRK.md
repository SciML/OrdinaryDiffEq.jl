```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqPDIRK

PDIRK methods are parallel DIRK methods.
SDIRK methods, or singly-diagonally implicit methods,
have to build and solve a factorize a Jacobian of the form `W = I-gammaJ` where `gamma` is dependent on the chosen method.
PDIRK methods use multiple different choices of `gamma`, i.e. `W_i = I-gamma_iJ`,
which are all used in the update process.
There are some advantages to this,
as no SDIRK method can be a higher order than 5,
while DIRK methods generally can have arbitrarily high order and lower error coefficients,
leading to lower errors at larger dt sizes.
With the right construction of the tableau,
these matrices can be factorized and the underlying steps can be computed in parallel,
which is why these are the parallel DIRK methods.

!!! warning "Experimental"

    `OrdinaryDiffEqPDIRK` is experimental,
    as there are no parallel DIRK tableaus that achieve good performance in the literature.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqPDIRK", "PDIRK44")
```

## Full list of solvers

```@docs
PDIRK44
```
