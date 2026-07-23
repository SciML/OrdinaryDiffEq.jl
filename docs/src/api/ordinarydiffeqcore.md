# OrdinaryDiffEqCore API

This page lists user-facing OrdinaryDiffEqCore API. The controller API has its
own page, and solver-author hooks are documented separately in the
[developer extension API](@ref Developer-Extension-API).

## Integrator objects

```@docs
OrdinaryDiffEqCore.ODEIntegrator
```

## Threading options

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
Lower-level nonlinear-solve hooks intended for solver authors are listed on the
[developer extension API](@ref Developer-Extension-API) page.

```@docs
OrdinaryDiffEqNonlinearSolve.NLAnderson
OrdinaryDiffEqNonlinearSolve.NLFunctional
OrdinaryDiffEqNonlinearSolve.NLNewton
OrdinaryDiffEqNonlinearSolve.HomotopyNonlinearSolveAlg
```
