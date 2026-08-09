# OrdinaryDiffEqCore API

This page lists user-facing OrdinaryDiffEqCore API. The controller API has its
own page, and solver-author hooks are documented separately in the
[developer extension API](https://docs.sciml.ai/OrdinaryDiffEq/stable/devtools/internals/public_api/).

## Integrator objects

```@docs
OrdinaryDiffEqCore.ODEIntegrator
```

## Threading options

These are the values accepted by the `thread` / `threading` keyword of the solvers
that support it. They are public API but are not exported by OrdinaryDiffEqCore or by
any umbrella package, so they must be written qualified as
`OrdinaryDiffEqCore.PolyesterThreads()` (after `using OrdinaryDiffEqCore`) or brought
into scope with `using OrdinaryDiffEqCore: BaseThreads, PolyesterThreads`. Naming them
bare, as v6 code did, gives an `UndefVarError`. `PolyesterThreads` additionally
requires `using Polyester`, which became a weak dependency in v7.

```@docs
OrdinaryDiffEqCore.AbstractThreadingOption
OrdinaryDiffEqCore.Sequential
OrdinaryDiffEqCore.BaseThreads
OrdinaryDiffEqCore.PolyesterThreads
OrdinaryDiffEqCore.isthreaded
```

## Automatic algorithm switching

```@docs
OrdinaryDiffEqCore.AutoAlgSwitch
OrdinaryDiffEqCore.AutoSwitch
```

## SSP helpers

```@docs
OrdinaryDiffEqCore.ssp_coefficient
```

## Limiter support

```@docs
OrdinaryDiffEqCore.has_stage_limiter
```

## Implicit method predictors

```@docs
OrdinaryDiffEqCore.Predictor
```

## Nonlinear solver algorithms

These OrdinaryDiffEqNonlinearSolve algorithms are public user API because they
are passed directly through implicit solver constructors as `nlsolve = ...`.

```@docs
OrdinaryDiffEqNonlinearSolve.NLAnderson
OrdinaryDiffEqNonlinearSolve.NLFunctional
OrdinaryDiffEqNonlinearSolve.NLNewton
OrdinaryDiffEqNonlinearSolve.NonlinearSolveAlg
OrdinaryDiffEqNonlinearSolve.HomotopyNonlinearSolveAlg
```
