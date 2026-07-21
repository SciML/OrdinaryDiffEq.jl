# Common Interface API

OrdinaryDiffEq re-exports the common problem, callback, solve, and automatic
differentiation interfaces needed to construct and solve ordinary differential
equations.

## Solve interface

```@docs
OrdinaryDiffEq.solve
OrdinaryDiffEq.solve!
OrdinaryDiffEq.init
OrdinaryDiffEq.step!
```

## Problem types

```@docs
OrdinaryDiffEq.ODEProblem
OrdinaryDiffEq.ODEFunction
OrdinaryDiffEq.ODESolution
OrdinaryDiffEq.SplitODEProblem
OrdinaryDiffEq.SplitFunction
OrdinaryDiffEq.SecondOrderODEProblem
OrdinaryDiffEq.DynamicalODEProblem
OrdinaryDiffEq.DAEProblem
OrdinaryDiffEq.DAEFunction
OrdinaryDiffEq.DAESolution
OrdinaryDiffEq.EnsembleProblem
```

### Ensemble context

```@docs; canonical=false
SciMLBase.EnsembleContext
```

## Callbacks

```@docs
OrdinaryDiffEq.CallbackSet
OrdinaryDiffEq.ContinuousCallback
OrdinaryDiffEq.DiscreteCallback
OrdinaryDiffEq.VectorContinuousCallback
```

## Solution and integrator utilities

```@docs
OrdinaryDiffEq.ReturnCode
OrdinaryDiffEq.ODEAliasSpecifier
OrdinaryDiffEq.add_tstop!
OrdinaryDiffEq.derivative_discontinuity!
OrdinaryDiffEq.reinit!
OrdinaryDiffEq.remake
OrdinaryDiffEq.set_proposed_dt!
OrdinaryDiffEq.successful_retcode
```

## Automatic differentiation

```@docs
OrdinaryDiffEq.AutoFiniteDiff
OrdinaryDiffEq.AutoForwardDiff
OrdinaryDiffEq.AutoSparse
```

## Default algorithm

```@docs
OrdinaryDiffEq.DefaultODEAlgorithm
```

## Interface modules

```@docs
OrdinaryDiffEq.SciMLBase
OrdinaryDiffEq.SciMLLogging
```
