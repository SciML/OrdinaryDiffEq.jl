# DiffEqBase API

This page lists the user-facing DiffEqBase API documented with OrdinaryDiffEq.
Solver-author hooks, callback machinery, and cache types are documented separately
in the [developer extension API](https://docs.sciml.ai/OrdinaryDiffEq/stable/devtools/internals/public_api/).

## Specialization levels

DiffEqBase implements parameter-container specialization levels owned and documented by
[SciMLBase](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/#specialization_levels).
They are re-exported from OrdinaryDiffEq for use with a plain `using OrdinaryDiffEq`.

`AutoDespecialize` accepts arbitrary parameter objects. During solver concretization,
DiffEqBase stores `p` in a `SciMLBase.DespecializedParameters` container with a stable
outer type. SciML function calls recover the original concrete parameter at a dynamic
function barrier, allowing precompiled solver code to be shared across parameter layouts.
`AutoSpecialize` retains its existing behavior and does not select this container.

`AutoRespecialize` is the constrained, non-dynamic policy formerly named
`AutoDePSpecialize`. Supported paths pack compatible parameters into an opaque container
and recover the original concrete type without dynamic dispatch. The deprecated
`AutoDePSpecialize` name remains available as an alias.

OrdinaryDiffEq re-exports `AutoDespecialize`, `AutoRespecialize`, and the deprecated
`AutoDePSpecialize` alias. Their canonical API documentation is on the linked SciMLBase
specialization-level page.

```@docs
OrdinaryDiffEq.AutoDespecialize
OrdinaryDiffEq.AutoRespecialize
OrdinaryDiffEq.AutoDePSpecialize
SciMLBase.AutoRespecialize
```

## Default callback behavior

```@docs
DiffEqBase.ODE_DEFAULT_ISOUTOFDOMAIN
DiffEqBase.ODE_DEFAULT_NORM
DiffEqBase.ODE_DEFAULT_PROG_MESSAGE
DiffEqBase.ODE_DEFAULT_UNSTABLE_CHECK
DiffEqBase.NAN_CHECK
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
