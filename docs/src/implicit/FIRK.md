```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqFIRK

Fully Implicit Runge-Kutta (FIRK) methods for stiff differential equations requiring very high accuracy. These methods solve a fully coupled implicit system at each timestep, providing superior accuracy and stability compared to diagonally implicit methods.

!!! warning "Real Numbers Only"
    
    FIRK methods should only be used for problems defined on real numbers, not complex numbers.

## Key Properties

FIRK methods provide:

  - **Highest-order implicit methods** (excluding extrapolation)
  - **Superior accuracy** for very low tolerance requirements (≤ 1e-9)
  - **A-stable and L-stable** behavior for stiff problems
  - **Higher order per stage** than SDIRK methods (order 2s+1 for s stages)
  - **Special geometric properties** (some methods are symplectic)
  - **Excellent for small to medium systems** with high accuracy requirements

## When to Use FIRK Methods

These methods are recommended for:

  - **Very low tolerance problems** (1e-9 and below) where accuracy is paramount
  - **Small to medium stiff systems** (< 200 equations)
  - **Problems requiring highest possible accuracy** for implicit methods
  - **Stiff problems** where SDIRK order limitations (max order 5) are insufficient
  - **Applications where computational cost is acceptable** for maximum accuracy

## Mathematical Background

RadauIIA methods are based on Gaussian collocation and achieve order 2s+1 for s stages, making them among the highest-order implicit methods available. They represent the ODE analog of Gaussian quadrature. For more details on recent advances in FIRK methods, see our paper: [High-Order Adaptive Time Stepping for the Incompressible Navier-Stokes Equations](https://arxiv.org/abs/2412.14362).

## Computational Considerations

### Advantages

  - **Higher accuracy per stage** than diagonal methods
  - **Better multithreading** for small systems due to larger linear algebra operations
  - **No order restrictions** like SDIRK methods (which max out at order 5)

### Disadvantages

  - **Limited to real-valued problems** - cannot be used for complex number systems
  - **Higher implementation complexity** compared to SDIRK methods

## Solver Selection Guide

### High accuracy requirements

  - **`AdaptiveRadau`**: **Recommended** - adaptive order method that automatically selects optimal order
  - **`RadauIIA5`**: 5th-order method, good balance of accuracy and efficiency
  - **`RadauIIA9`**: 9th-order method for extremely high accuracy requirements
  - **`RadauIIA3`**: 3rd-order method for moderate accuracy needs

### System size considerations

  - **Systems < 200**: FIRK methods are competitive due to better multithreading
  - **Systems > 200**: Consider SDIRK or BDF methods instead

## Performance Guidelines

  - **Best for tolerances ≤ 1e-9** where high accuracy justifies the cost
  - **Most efficient on small to medium systems** where linear algebra cost is manageable
  - **Should be tested against** parallel implicit extrapolation methods which specialize in similar regimes
  - **Compare with** high-order SDIRK methods for borderline cases

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqFIRK", "RadauIIA5")
```

## Full list of solvers

```@docs
RadauIIA3
RadauIIA5
RadauIIA9
```
