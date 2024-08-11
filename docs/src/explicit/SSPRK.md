# OrdinaryDiffEqSSPRK

Solvers specialized on certain PDE types.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqSSPRK", "SSPRK22")
```

## Full list of solvers

```@docs
SSPRK22
SSPRK33
SSPRK53
KYKSSPRK42
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
SHLDDRK_2N
SHLDDRK52
```