```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqPRK

Parallel Runge-Kutta (PRK) methods are explicit solvers specifically designed to exploit parallelism by making multiple independent evaluations of the ODE function `f` simultaneously. These methods are optimized for parallel computing environments where function evaluations can be distributed across multiple processors.

!!! warning "Research and Development"
    
    These methods are currently in research and development and not intended for general use.

## Key Properties

PRK methods provide:

  - **Explicit parallelism** in function evaluations within each timestep
  - **Fixed processor count optimization** for specific parallel architectures
  - **Independent stage evaluations** that can run simultaneously
  - **Maintained accuracy** while achieving parallel speedup
  - **Specialized tableaus** designed for parallel efficiency

## When to Use PRK Methods

These methods are recommended for:

  - **Parallel computing environments** with multiple processors available
  - **Expensive function evaluations** that benefit from parallelization
  - **Systems where function evaluation dominates** computational cost
  - **Applications with fixed parallel architecture** (e.g., exactly 2 or 5 processors)
  - **Problems where parallel speedup** outweighs method overhead

## Important Considerations

### Parallel Requirements

  - **Requires multiple processors** to achieve benefits
  - **Function evaluations must be parallelizable** (no data dependencies)
  - **Parallel overhead** must be less than speedup gains
  - **Fixed processor count** optimization may not match available hardware

### When NOT to Use

  - **Sequential computing** environments
  - **Cheap function evaluations** where parallel overhead dominates
  - **Memory-bound problems** where parallelism doesn't help
  - **Variable processor availability** scenarios
  - **Large systems** where LU factorization of implicit steps parallelizes efficiently (around 200×200 matrices and larger on modern processors)

## Mathematical Background

PRK methods rearrange traditional Runge-Kutta tableaus to allow stage evaluations to be computed independently and simultaneously. The specific processor count determines the tableau structure and achievable parallelism.

## Solver Selection Guide

### Available methods

  - **`KuttaPRK2p5`**: Fifth-order method optimized for 2 processors

### Usage considerations

  - **Best with exactly 2 processors** for KuttaPRK2p5
  - **Function evaluation must support** parallel execution
  - **Test parallel efficiency** against sequential high-order methods
  - **Consider problem-specific** parallel architecture

## Performance Guidelines

  - **Measure actual speedup** vs sequential methods on target hardware
  - **Account for parallel overhead** in performance comparisons
  - **Consider memory bandwidth** limitations in parallel environments
  - **Compare against** other parallelization strategies (e.g., spatial domain decomposition)

## Alternative Parallelization Approaches

For most problems, consider these alternatives:

  - **Spatial domain decomposition** for PDE problems
  - **Multiple trajectory parallelism** for Monte Carlo simulations
  - **Vectorized operations** within function evaluations
  - **High-order sequential methods** with better single-thread performance

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqPRK", "KuttaPRK2p5")
```

## Full list of solvers

```@docs
KuttaPRK2p5
```
