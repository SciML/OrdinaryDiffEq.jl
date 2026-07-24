# DiffEqBase API

This page lists the user-facing DiffEqBase API documented with OrdinaryDiffEq.
Solver-author hooks, callback machinery, and cache types are documented separately
in the [developer extension API](https://docs.sciml.ai/OrdinaryDiffEq/stable/devtools/internals/public_api/).

## Specialization levels

DiffEqBase implements the ODE-path wrapping for the `AutoDePSpecialize`
specialization level, which is owned and documented by
[SciMLBase](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/#specialization_levels).
It is re-exported from OrdinaryDiffEq, so
`ODEProblem{true, AutoDePSpecialize}(f!, u0, tspan, p)` works from a plain
`using OrdinaryDiffEq`. It packs an `isbits` parameter into a
`RespecializeParams.OpaqueParams` container so one precompiled solve is shared
across all parameter struct types; recover the original payload from
`sol.prob.p` with `RespecializeParams.unpack(sol.prob.p, typeof(p))`.

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
