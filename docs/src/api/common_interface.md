# Common Interface API

OrdinaryDiffEq re-exports the common problem, callback, solve, and automatic
differentiation interfaces needed to construct and solve ordinary differential
equations. The exact set of
re-exported SciMLBase names, and the rule that decides it, are given in
[Re-exported SciMLBase API](@ref).

## Solve interface

The generic solve interface is selected by the problem and algorithm types; it
does not require callers to depend on an OrdinaryDiffEq implementation module.
Use `solve` for a complete solution, `init` when an integrator must be inspected
or advanced incrementally, `step!` to advance one accepted step, and `solve!` to
finish an initialized integrator. Algorithm-specific constructors are the only
solver-family dependency in the normal user workflow.

```@autodocs
Modules = [CommonSolve, SciMLBase]
Public = true
Private = false
Filter = x -> x in [
    CommonSolve.solve,
    CommonSolve.solve!,
    CommonSolve.init,
    SciMLBase.step!
]
```

## Problem types

Construct a problem with the function signature required by its problem type:
out-of-place functions return the next state, while in-place functions write to
their first argument and return `nothing`. The problem owns the time span and
parameters; solver options belong in `solve` or `init` keyword arguments. Use
`remake` to derive a new problem while preserving the original function and
metadata.

```@autodocs
Modules = [SciMLBase]
Public = true
Private = false
Filter = x -> x in [
    SciMLBase.ODEProblem,
    SciMLBase.ODEFunction,
    SciMLBase.ODENLStepData,
    SciMLBase.ODESolution,
    SciMLBase.SplitODEProblem,
    SciMLBase.SplitFunction,
    SciMLBase.SecondOrderODEProblem,
    SciMLBase.DynamicalODEProblem,
    SciMLBase.DAEProblem,
    SciMLBase.DAEFunction,
    SciMLBase.DAESolution,
    SciMLBase.DiscreteProblem,
    SciMLBase.DiscreteFunction,
    SciMLBase.EnsembleProblem,
    SciMLBase.DynamicalODEFunction,
    SciMLBase.EnsembleContext
]
```

### Ensemble context

## Callbacks

Callbacks are supplied through the problem or solver call. A continuous
condition returns a signed value whose zero crossing triggers `affect!`; a
discrete condition returns a Boolean. Callback effects mutate only through the
integrator interface, and termination should use the generic `terminate!`
operation rather than a solver-specific field update.

```@autodocs
Modules = [SciMLBase]
Public = true
Private = false
Filter = x -> x in [
    SciMLBase.CallbackSet,
    SciMLBase.ContinuousCallback,
    SciMLBase.DiscreteCallback,
    SciMLBase.VectorContinuousCallback
]
```

## Solution and integrator utilities

Time stops must lie in the integration direction and are reached exactly by the
generic integrator. `successful_retcode` is the portable way to test completion;
applications should inspect `ReturnCode` rather than relying on a solver's
internal status representation.

```@autodocs
Modules = [SciMLBase]
Public = true
Private = false
Filter = x -> x in [
    SciMLBase.ReturnCode,
    SciMLBase.ODEAliasSpecifier,
    SciMLBase.add_tstop!,
    SciMLBase.derivative_discontinuity!,
    SciMLBase.reinit!,
    SciMLBase.remake,
    SciMLBase.set_proposed_dt!,
    SciMLBase.successful_retcode,
    SciMLBase.has_global_error,
    SciMLBase.add_saveat!,
    SciMLBase.auto_dt_reset!
]
```

## Automatic differentiation

```@autodocs
Modules = [ADTypes]
Public = true
Private = false
Filter = x -> nameof(x) in (
    :AutoFiniteDiff,
    :AutoForwardDiff,
    :AutoSparse,
    :AbstractADType,
    :AbstractSparsityDetector,
    :AbstractColoringAlgorithm,
    :mode,
    :AbstractMode,
    :jacobian_sparsity,
    :hessian_sparsity,
    :KnownJacobianSparsityDetector,
    :KnownHessianSparsityDetector,
    :ForwardMode,
    :ReverseMode,
    :ForwardOrReverseMode,
    :SymbolicMode
)
```

The solver-side differentiation interface consumes ADTypes backends and
SciMLOperators for matrix-free Jacobian and linear-solver integrations.

```@autodocs
Modules = [SciMLOperators]
Public = true
Private = false
Filter = x -> x === SciMLOperators.AbstractSciMLOperator
```

## Default algorithm

```@docs
OrdinaryDiffEqDefault.DefaultODEAlgorithm
OrdinaryDiffEqDefault.DefaultImplicitODEAlgorithm
```
