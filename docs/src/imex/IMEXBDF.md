```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqBDF

IMEX BDF (Implicit-Explicit Backward Differentiation Formula) methods for stiff differential equations that can be split into stiff and non-stiff components. These methods apply implicit BDF schemes to the stiff part while treating the non-stiff part explicitly, providing efficient handling of problems with mixed stiffness characteristics.

## Key Properties

IMEX BDF methods provide:

  - **Implicit-explicit splitting** for mixed stiffness problems
  - **BDF stability** for the stiff component with A-stable and L-stable behavior
  - **Explicit treatment** of non-stiff terms avoiding unnecessary computational cost
  - **High-order accuracy** up to 4th order for both components
  - **Efficient for large systems** where full implicit treatment is expensive
  - **Natural for operator splitting** problems

## When to Use IMEX BDF Methods

These methods are recommended for:

  - **Reaction-diffusion systems** where reaction terms are stiff and diffusion is moderate
  - **Convection-diffusion problems** with stiff source terms and explicit convection
  - **Parabolic PDEs** where diffusion operators are naturally split from other terms
  - **Problems with natural stiffness separation** where some terms require implicit treatment
  - **Large-scale systems** where full implicit methods are computationally prohibitive
  - **Applications requiring operator splitting** methodology

## Mathematical Background

IMEX BDF methods split the ODE system `du/dt = f(u,t)` into:
`du/dt = f₁(u,t) + f₂(u,t)`

where:

  - `f₁(u,t)` contains stiff terms (treated implicitly with BDF)
  - `f₂(u,t)` contains non-stiff terms (treated explicitly)

This splitting must be chosen carefully to ensure both stability and efficiency.

## Problem Splitting Requirements

These methods require a `SplitODEProblem` formulation where:

  - **First function** `f₁` should contain stiff, implicit terms
  - **Second function** `f₂` should contain non-stiff, explicit terms
  - **Splitting strategy** significantly affects method performance
  - **Stiffness characteristics** should align with implicit/explicit treatment

## Solver Selection Guide

### IMEX Multistep Methods

  - **`SBDF2`**: **Recommended** - Second-order IMEX BDF method, good balance of accuracy and stability
  - **`SBDF3`**: Third-order method for higher accuracy requirements
  - **`SBDF4`**: Fourth-order method for maximum accuracy in IMEX BDF family
  - **`SBDF`**: Adaptive order method (experimental)

### IMEX SDIRK Methods

  - **`IMEXEuler`**: First-order method for simple problems or debugging
  - **`IMEXEulerARK`**: Alternative first-order formulation

## Performance Guidelines

### When IMEX BDF methods excel

  - **Natural stiffness separation** where splitting is obvious
  - **Large systems** where full implicit treatment is expensive
  - **Parabolic PDEs** with natural operator splitting
  - **Reaction-diffusion problems** with well-separated timescales
  - **Problems where implicit component** has efficient linear algebra

### Splitting strategy considerations

  - **Identify stiff vs non-stiff terms** based on eigenvalue analysis
  - **Linear stiff terms** work well in implicit component
  - **Nonlinear non-stiff terms** are suitable for explicit treatment
  - **Test different splittings** to optimize performance

## Alternative Approaches

Consider these alternatives:

  - **Full implicit methods** (BDF, SDIRK) if splitting is unclear or ineffective
  - **Standard IMEX Runge-Kutta** methods for different accuracy/efficiency trade-offs
  - **Exponential integrators** for linear stiff problems with nonlinear non-stiff terms
  - **Rosenbrock methods** for moderately stiff problems without natural splitting

## Usage Considerations

  - **Careful splitting design** is crucial for method effectiveness
  - **Stability analysis** should verify that explicit treatment doesn't introduce instabilities
  - **Timestep restrictions** may apply to the explicit component
  - **Linear algebra efficiency** in the implicit component affects overall performance

```@eval
imex_first_steps = evalfile("./common_imex_first_steps.jl")
imex_first_steps("OrdinaryDiffEqBDF", "SBDF2")
```

## Full list of solvers

### IMEX Multistep

```@docs
SBDF
SBDF2
SBDF3
SBDF4
```

### IMEX SDIRK

Note that Implicit Euler is the 1st order BDF method, and is thus implemented here using
the same machinery.

```@docs
IMEXEuler
IMEXEulerARK
```
