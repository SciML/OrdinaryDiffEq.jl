```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqSSPRK

SSPRK methods are Runge-Kutta methods which support the "strongly preserving property" (SSP).
They are designed for the use in discretizations of hyperbolic partial differential equations and conservation laws
as they have extra stability properties ( e.g., stability with respect to total variation, the maximum norm, or other convex functionals)
when step-size restrictions are respected.
In particular, these properties are granted if the step-size is kept to a level where the CFL coefficients are less than the SSP coefficient.

Note that for SSPRK methods, a algorithm utility `OrdinaryDiffEqCore.ssp_coefficient(alg)` exists that allows for querying the SSP coefficient for use in step size calculations.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqSSPRK", "SSPRK432")
```

## Full list of solvers

```@docs
SSPRK22
SSPRK33
SSPRK53
KYKSSPRK42
KYK2014DGSSPRK_3S2
SSPRK53_2N1
SSPRK53_2N2
SSPRK53_H
SSPRK63
SSPRK73
SSPRK83
SSPRK43
SSPRK432
SSPRKMSVS43
SSPRKMSVS32
SSPRK932
SSPRK54
SSPRK104
```
