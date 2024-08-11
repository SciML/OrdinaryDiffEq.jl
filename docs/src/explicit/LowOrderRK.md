# OrdinaryDiffEqLowOrderRK

If [`OrdinaryDiffEqTsit5`](@ref OrdinaryDiffEqTsit5) is not working well for your non-stiff problem at default and higher tolerance,
it can be worthwhile to explore the options in this package.
In particular, when more robust error control is required, [`BS5`](@ref) is a good choice.
If at moderate tolerances and the interpolation error is very important,
consider the [`OwrenZen5`](@ref) method.
For fast solving at higher tolerances, we recommend [`BS3`](@ref),
or [`OwrenZen3`](@ref)if the interpolation error is important.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqLowOrderRK", "BS3")
```

## Full list of solvers

```@docs
Euler
Heun
Ralston
Midpoint
RK4
BS3
OwrenZen3
OwrenZen4
OwrenZen5
BS5
DP5
Anas5
RKO65
FRK65
RKM
MSRK5
MSRK6
PSRK4p7q6
PSRK3p5q4
PSRK3p6q5
Stepanov5
SIR54
Alshina2
Alshina3
Alshina6
```
