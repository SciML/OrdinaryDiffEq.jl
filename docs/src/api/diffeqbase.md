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
```

```@autodocs
Modules = [SciMLBase]
Public = true
Private = false
Filter = x -> x === SciMLBase.AutoRespecialize
```

### What each level reuses

For an in-place `ODEProblem` solved with OrdinaryDiffEq, the table shows which changes
between two problems keep the integrator type, so the second solve reuses the solver
compiled for the first. "Recompiles" means the second problem compiles a new solver.

| Level | New `f` | New parameter type | New callback type | New `sys` type |
|---|---|---|---|---|
| `FullSpecialize` | recompiles | recompiles | recompiles | recompiles |
| `AutoSpecialize` | reused | recompiles | reused on Julia 1.12+ | recompiles |
| `AutoDespecialize` | reused | reused | reused on Julia 1.12+ | recompiles |
| `AutoRespecialize` | reused | reused among `isbits` types and among other types | recompiles | recompiles |
| `FunctionWrapperSpecialize` | reused | recompiles | recompiles | recompiles |
| `NoSpecialize` | reused | recompiles | reused on Julia 1.12+ | recompiles |

The automatic levels only wrap `f` when `u0` is an array of unitless numbers, not a
`SubArray`, `t` is unitless and `f` is not a SciMLOperator. With solvers that do not use
ForwardDiff, such as `Tsit5`, `f` must also have no analytic `jac`. Other problems
specialize on `f` as `FullSpecialize` does, and under `AutoRespecialize` on the parameter
type as well.

`AutoRespecialize` packs `isbits` parameters into one container type and all other
parameters into another, so a change from one kind to the other compiles a new solver. It
leaves the parameters alone when `f` has a `sys`, where every new parameter type
recompiles.

Callback types are erased only on Julia 1.12 and later. On older versions every level
compiles a new solver for each new callback type. Where callbacks are erased, the first
solve with a callback still compiles the callback handling once, and later callback types
compile only their own condition and affect functions.

A new `sys` type, the symbolic container used for symbolic indexing, compiles a new solver
at every level.

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
