```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqSSPRK

Strong Stability Preserving Runge-Kutta (SSPRK) methods are specialized explicit Runge-Kutta methods designed to preserve important stability properties of the underlying spatial discretization when applied to hyperbolic partial differential equations and conservation laws.

## Key Properties

SSPRK methods provide:

  - **Strong stability preservation** for convex functionals (total variation, maximum norm, entropy)
  - **Optimal SSP coefficients** allowing larger stable timesteps
  - **Non-oscillatory behavior** crucial for hyperbolic PDEs and conservation laws
  - **High-order accuracy** while maintaining monotonicity properties
  - **Specialized variants** for different orders and storage requirements

## When to Use SSPRK Methods

SSPRK methods are essential for:

  - **Hyperbolic partial differential equations** (Euler equations, shallow water, etc.)
  - **Conservation laws** where preserving physical bounds is critical
  - **Discontinuous Galerkin methods** and other high-order spatial discretizations
  - **Problems requiring monotonicity preservation** or total variation stability
  - **Shock-capturing schemes** where spurious oscillations must be avoided
  - **Astrophysical simulations** and computational fluid dynamics

## SSP Coefficient and CFL Conditions

The SSP coefficient determines the maximum allowable timestep for stability preservation. Use `OrdinaryDiffEqCore.ssp_coefficient(alg)` to query this value for step size calculations. The timestep must satisfy `dt ≤ CFL * dx / max_wavespeed` where CFL ≤ SSP coefficient.

## Solver Selection Guide

### Second-order methods

  - **`SSPRK22`**: Two-stage, second-order (SSP coefficient = 1)
  - **`SSPRKMSVS32`**: Three-step multistep variant

### Third-order methods

  - **`SSPRK33`**: Three-stage, third-order, optimal (SSP coefficient = 1)
  - **`SSPRK53`**: Five-stage, third-order, higher SSP coefficient
  - **`SSPRK63`**, **`SSPRK73`**, **`SSPRK83`**: More stages for larger SSP coefficients
  - **`SSPRK43`**: Four-stage with embedded error estimation
  - **`SSPRK432`**: Low-storage variant

### Fourth-order methods

  - **`SSPRK54`**: Five-stage, fourth-order
  - **`SSPRK104`**: Ten-stage, fourth-order, large SSP coefficient

### Low-storage variants

  - **`SSPRK53_2N1`**, **`SSPRK53_2N2`**, **`SSPRK53_H`**: Two-register storage schemes

### Discontinuous Galerkin optimized

  - **`KYKSSPRK42`**: Optimized for DG spatial discretizations
  - **`KYK2014DGSSPRK_3S2`**: Specialized DG method

### Adaptive methods

  - **`SSPRK432`**: Third-order with error control
  - **`SSPRK932`**: High-stage adaptive method
  - **`SSPRKMSVS43`**: Multistep adaptive variant

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqSSPRK", "SSPRK22")
```

## Full list of solvers

```@docs
SSPRK22
SSPRK33
SSPRK53
KYKSSPRK42
KYK2014DGSSPRK_3S2
SSPRK53_2N1
SSPRK53_2N2
SSPRK53_H
SSPRK63
SSPRK73
SSPRK83
SSPRK43
SSPRK432
SSPRKMSVS43
SSPRKMSVS32
SSPRK932
SSPRK54
SSPRK104
```
