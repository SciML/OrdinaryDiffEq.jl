# DiffEqBase API

This page lists the user-facing DiffEqBase API documented with OrdinaryDiffEq.
Solver-author hooks, callback machinery, and cache types are documented separately
in the [developer extension API](@ref Developer-Extension-API).

## Default callback behavior

```@docs
DiffEqBase.ODE_DEFAULT_ISOUTOFDOMAIN
DiffEqBase.ODE_DEFAULT_NORM
DiffEqBase.ODE_DEFAULT_PROG_MESSAGE
DiffEqBase.ODE_DEFAULT_UNSTABLE_CHECK
```

## Runge-Kutta tableau types

```@docs
DiffEqBase.Tableau
DiffEqBase.ODERKTableau
DiffEqBase.ExplicitRKTableau
DiffEqBase.ImplicitRKTableau
```

## Cost and convergence helpers

```@docs
DiffEqBase.ConvergenceSetup
DiffEqBase.DECostFunction
```

## DAE initialization

```@docs
DiffEqBase.DefaultInit
DiffEqBase.BrownFullBasicInit
DiffEqBase.ShampineCollocationInit
```

## Sensitivity passthrough

```@docs
DiffEqBase.SensitivityADPassThrough
```
