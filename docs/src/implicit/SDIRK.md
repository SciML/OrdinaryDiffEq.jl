```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqSDIRK

Singly Diagonally Implicit Runge-Kutta (SDIRK) methods are a family of implicit Runge-Kutta methods designed for solving stiff ordinary differential equations. These methods are particularly effective for problems where explicit methods become unstable due to stiffness.

## Key Properties

SDIRK methods have several important characteristics:

  - **A-stable and L-stable**: Can handle highly stiff problems without numerical instability
  - **Stiffly accurate**: Many SDIRK methods provide additional numerical stability for stiff problems
  - **Diagonally implicit structure**: The implicit system only requires solving a sequence of nonlinear equations rather than a large coupled system
  - **Good for moderate to large systems**: More efficient than fully implicit RK methods for many problems

## When to Use SDIRK Methods

SDIRK methods are recommended for:

  - **Stiff differential equations** where explicit methods fail or require very small timesteps
  - **Problems requiring good stability properties** at moderate to high tolerances
  - **Systems where Rosenbrock methods** (which require Jacobians) are not suitable or available
  - **IMEX problems** using the KenCarp family, which can split stiff and non-stiff terms

## Solver Selection Guide

### High tolerance (>1e-2)

  - **`TRBDF2`**: Second-order A-B-L-S-stable method, good for oscillatory problems

### Medium tolerance (1e-8 to 1e-2)

  - **`KenCarp4`**: Fourth-order method with excellent stability, good all-around choice
  - **`KenCarp47`**: Seventh-stage fourth-order method, enhanced stability
  - **`Kvaerno4`** or **`Kvaerno5`**: High-order stiffly accurate methods

### Low tolerance (<1e-8)

  - **`Kvaerno5`**: Fifth-order stiffly accurate method for high accuracy
  - **`KenCarp5`**: Fifth-order method with splitting capabilities

### Special Cases

  - **`ImplicitEuler`**: First-order method, only recommended for problems with discontinuities or when `f` is not differentiable
  - **`Trapezoid`**: Second-order symmetric method, reversible but not symplectic. Good for eliminating damping often seen with L-stable methods
  - **`ImplicitMidpoint`**: Second-order A-stable symplectic method for energy-preserving systems
  - **`SSPSDIRK2`**: Strong stability preserving variant for problems requiring monotonicity preservation

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqSDIRK", "KenCarp4")
```

## Full list of solvers

```@docs
ImplicitEuler
ImplicitMidpoint
Trapezoid
TRBDF2
SDIRK2
SDIRK22
SSPSDIRK2
Kvaerno3
KenCarp3
CFNLIRK3
Cash4
SFSDIRK4
SFSDIRK5
SFSDIRK6
SFSDIRK7
SFSDIRK8
Hairer4
Hairer42
Kvaerno4
Kvaerno5
```

### IMEX SDIRK

These methods support `SplitODEProblem` for implicit-explicit (IMEX) integration,
where the stiff part is treated implicitly and the non-stiff part is treated explicitly.

```@docs
KenCarp4
KenCarp47
KenCarp5
KenCarp58
```

### Higher-Order ESDIRK

Higher-order ESDIRK (Explicit first stage Singly Diagonally Implicit Runge-Kutta) methods.
These are high-order L-stable implicit methods where the first stage is explicit.

!!! note
    These methods do not support `SplitODEProblem`. For IMEX integration with split
    problems, use the KenCarp methods above instead.

```@docs
ESDIRK54I8L2SA
ESDIRK436L2SA2
ESDIRK437L2SA
ESDIRK547L2SA2
ESDIRK659L2SA
```
