# OrdinaryDiffEqCore API

This page lists user-facing OrdinaryDiffEqCore API. The controller API has its
own page, and solver-author hooks are documented separately in the
[developer extension API](https://docs.sciml.ai/OrdinaryDiffEq/stable/devtools/internals/public_api/).

## Integrator objects

```@docs
OrdinaryDiffEqCore.ODEIntegrator
```

## Threading options

These are the values accepted by the `threading` keyword of the solvers that expose
independent internal work — the extrapolation methods, the parallel Runge-Kutta
methods (`KuttaPRK2p5`, `PDIRK44`) and the parallel-stage FIRK methods:

```julia
using OrdinaryDiffEqExtrapolation

prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
solve(
    prob,
    ExtrapolationMidpointDeuflhard(;
        threading = OrdinaryDiffEqExtrapolation.BaseThreads()
    )
)
```

`Sequential`, `BaseThreads` and `PolyesterThreads` are public API of
OrdinaryDiffEqCore and are declared public by OrdinaryDiffEq and by each
sublibrary that takes `threading`, so they are reachable qualified through
whichever of those you already have loaded. They are deliberately not exported —
write them qualified, or bring them into scope with
`using OrdinaryDiffEqCore: BaseThreads, PolyesterThreads`. `PolyesterThreads`
additionally requires `using Polyester`, which became a weak dependency in v7.

The `thread` keyword of the FastBroadcast-based solvers is a different option that
takes `FastBroadcast.Serial()` or `FastBroadcast.Threaded()`; these types are not
interchangeable with the ones below.

```@docs
OrdinaryDiffEqCore.AbstractThreadingOption
OrdinaryDiffEqCore.Sequential
OrdinaryDiffEqCore.BaseThreads
OrdinaryDiffEqCore.PolyesterThreads
OrdinaryDiffEqCore.isthreaded
```

The umbrella package exposes the same documented threading options as qualified
public names:

```@docs
OrdinaryDiffEq.Sequential
OrdinaryDiffEq.BaseThreads
OrdinaryDiffEq.PolyesterThreads
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
