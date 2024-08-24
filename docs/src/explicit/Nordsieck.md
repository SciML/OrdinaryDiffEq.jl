```@meta
CollapsedDocStrings = true
```
# OrdinaryDiffEqNordsieck

The Nordsieck form is an alternative representation of multistep methods which,
instead of representing and saving past step values in a history vector,
it uses a derivative list (like a Taylor expansion) for the computation of the next point.
The Nordsieck form was pioneered by early implementations of BDF methods such LSODE, VODE, and finally CVODE.
It can have some advantages in terms of restartability as the full Nordsieck vector can be instantiated given only the information of f and its derivatives after discontinuities,
but the higher derivative representations can also introduce numerical instabilities of their own.

The Nordsieck implementations here are considered experimental implementations of the LSODE non-fixed leading coefficient form
and are generally considered inferior to the fixed-leading history-based BDF implementation of FBDF, and thus for all standard usage we recommend FBDF.
However, this algorithm is kept for experimental research and development purposes with the possibility of one day becoming a more discontinuity-aware BDF implementation.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqNordsieck", "AN5")
```

## Full list of solvers

```@docs
AN5
JVODE
JVODE_Adams
JVODE_BDF
```
