```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqStabilizedIRK

Stabilized Implicit Runge-Kutta (IMEX) methods combine the benefits of stabilized explicit methods with implicit treatment of problematic eigenvalues. These IMEX schemes are designed for problems where the Jacobian has both large real eigenvalues (suitable for explicit stabilized methods) and large complex eigenvalues (requiring implicit treatment).

## Key Properties

Stabilized IRK methods provide:

  - **IMEX formulation** treating different stiffness components appropriately
  - **Large stability regions** for real eigenvalues via explicit stabilized schemes
  - **Implicit treatment** of complex eigenvalues for unconditional stability
  - **Efficient handling** of mixed stiffness characteristics
  - **Splitting-based approach** requiring `SplitODEProblem` formulation

## When to Use Stabilized IRK Methods

These methods are recommended for:

  - **Mixed stiffness problems** with both real and complex eigenvalues
  - **Parabolic PDEs with convection** where diffusion and advection have different scales
  - **Reaction-diffusion systems** with stiff reactions and moderate diffusion
  - **Problems where pure explicit stabilized methods fail** due to complex eigenvalues
  - **Large-scale systems** where full implicit methods are too expensive

## Mathematical Background

Standard stabilized explicit methods (like RKC, ROCK) achieve large stability regions along the negative real axis but struggle with complex eigenvalues. Stabilized IRK methods address this by:

 1. **Explicit stabilized treatment** for large real eigenvalues
 2. **Implicit treatment** for complex eigenvalues
 3. **IMEX coupling** to maintain overall stability and accuracy

## Problem Splitting Requirements

These methods require a `SplitODEProblem` where:

  - **First component** contains terms with large real eigenvalues (explicit treatment)
  - **Second component** contains terms with complex eigenvalues (implicit treatment)
  - **Splitting design** is crucial for method performance

## Spectral Radius Estimation

Users can supply an upper bound on the spectral radius:

```julia
eigen_est = (integrator) -> integrator.eigen_est = upper_bound
```

This bound applies to the explicit component of the split problem.

## Solver Selection Guide

### Available methods

  - **`IRKC`**: Implicit Runge-Kutta-Chebyshev method for mixed stiffness problems

### Usage considerations

  - **Requires careful splitting** of the problem components
  - **Spectral radius estimation** needed for explicit component
  - **Test splitting strategies** for optimal performance
  - **Compare with pure implicit** or explicit stabilized alternatives

## Performance Guidelines

### When IMEX stabilized methods excel

  - **Mixed eigenvalue distribution** (both real and complex)
  - **Moderate to large systems** where splitting is natural
  - **Problems where neither pure explicit nor implicit** methods are ideal

### Splitting strategy considerations

  - **Identify dominant eigenvalue types** in different terms
  - **Real-dominated terms** → explicit component
  - **Complex-dominated terms** → implicit component
  - **Test different splittings** for best performance

## Alternative Approaches

Consider these alternatives:

  - **Pure implicit methods** (BDF, SDIRK, Rosenbrock) for highly stiff problems
  - **Explicit stabilized methods** (ROCK, RKC) if complex eigenvalues are small
  - **Standard IMEX methods** for natural explicit/implicit splitting

```@eval
imex_first_steps = evalfile("./common_imex_first_steps.jl")
imex_first_steps("OrdinaryDiffEqStabilizedIRK", "IRKC")
```

## Full list of solvers

```@docs
IRKC
```
