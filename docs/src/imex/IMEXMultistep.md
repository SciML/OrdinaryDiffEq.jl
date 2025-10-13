```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqIMEXMultistep

Standard low-order IMEX (Implicit-Explicit) multistep methods for problems that can be split into stiff and non-stiff components. These are widely used classical methods in partial differential equation applications, providing simple and reliable IMEX integration with moderate accuracy requirements.

## Key Properties

IMEX Multistep methods provide:

  - **Standard IMEX formulations** commonly used in PDE applications
  - **Low-order accuracy** (typically 2nd order) with good stability
  - **Simple implementation** and well-understood behavior
  - **Explicit treatment** of non-stiff terms with implicit handling of stiff components
  - **Fixed timestep requirements** due to multistep nature
  - **Efficient for large-scale problems** where splitting is natural

## When to Use IMEX Multistep Methods

These methods are recommended for:

  - **Classical PDE applications** where standard IMEX methods are established
  - **Reaction-diffusion systems** with natural explicit/implicit splitting
  - **Convection-diffusion problems** where convection is explicit and diffusion implicit
  - **Large-scale spatial discretizations** where simple, efficient methods are preferred
  - **Applications prioritizing robustness** over high-order accuracy
  - **Problems with natural operator splitting** methodology

## Mathematical Background

IMEX multistep methods treat the split system:
`du/dt = f₁(u,t) + f₂(u,t)`

using:

  - **Implicit multistep schemes** (like Crank-Nicolson) for stiff terms f₁
  - **Explicit multistep schemes** (like Adams-Bashforth) for non-stiff terms f₂

This combination provides stability for stiff components while maintaining efficiency for non-stiff parts.

## Problem Splitting Requirements

These methods require a `SplitODEProblem` where:

  - **First function** `f₁` contains stiff terms requiring implicit treatment
  - **Second function** `f₂` contains non-stiff terms suitable for explicit treatment
  - **Splitting should align** with the natural time scale separation
  - **Linear stiff terms** work particularly well with these methods

## Solver Selection Guide

### Available Methods

  - **`CNAB2`**: **Recommended** - Crank-Nicolson Adams-Bashforth 2nd order method
  - **`CNLF2`**: Crank-Nicolson Leap-Frog 2nd order method

### Method characteristics

  - **`CNAB2`**: Most commonly used, good stability and accuracy balance
  - **`CNLF2`**: Alternative formulation, may have different stability properties

## Performance Guidelines

### When IMEX Multistep methods excel

  - **PDE problems** with established IMEX splitting practices
  - **Large spatial discretizations** where method efficiency matters more than high accuracy
  - **Problems with linear stiff terms** that are efficiently handled implicitly
  - **Applications requiring consistent timesteps** (no adaptive timestepping)
  - **Well-conditioned problems** where simple methods suffice

### Splitting strategy considerations

  - **Linear diffusion terms** → implicit component (f₁)
  - **Nonlinear convection/reaction** → explicit component (f₂) if not too stiff
  - **Source terms** → choose based on stiffness characteristics
  - **Boundary conditions** → often naturally handled in implicit component

## Limitations and Considerations

### Method limitations

  - **Fixed timestep required** - no adaptive timestepping capabilities
  - **Low order only** - maximum 2nd order accuracy
  - **Startup procedures** needed for multistep methods
  - **Limited stability analysis** compared to modern IMEX-RK methods

### When to consider alternatives

  - **Higher accuracy needs**: Use IMEX-RK or higher-order IMEX methods
  - **Adaptive timestepping**: Use IMEX-RK or ARK methods
  - **Complex stability requirements**: Use more sophisticated IMEX schemes
  - **Very stiff problems**: Consider fully implicit methods

## Alternative Approaches

Consider these alternatives:

  - **IMEX Runge-Kutta** methods for adaptive timestepping and higher order
  - **IMEX BDF methods** for better stability properties and higher accuracy
  - **Fully implicit methods** if splitting is not beneficial
  - **Exponential integrators** for linear stiff problems

## Classical Applications

These methods are standard in:

  - **Computational fluid dynamics** for incompressible Navier-Stokes equations
  - **Atmospheric modeling** for advection-diffusion-reaction systems
  - **Ocean modeling** for transport equations with diffusion
  - **Astrophysical simulations** for multiphysics problems

```@eval
imex_first_steps = evalfile("./common_imex_first_steps.jl")
imex_first_steps("OrdinaryDiffEqIMEXMultistep", "CNAB2")
```

## Full list of solvers

```@docs
CNAB2
CNLF2
```
