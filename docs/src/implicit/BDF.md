```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqBDF

Backward Differentiation Formula (BDF) methods are multistep implicit methods specifically designed for solving large stiff systems of differential equations. They are the preferred choice for very large systems (>1000 equations) where other implicit methods become computationally expensive.

## Key Properties

BDF methods offer:

  - **Excellent efficiency for large systems** (>1000 ODEs)
  - **L-stable behavior** for orders 1 and 2 only
  - **Adaptive order and stepsize** control for optimal performance
  - **Alpha-stability** for higher orders (but less stable than L-stable methods for problems with large complex eigenvalues)

## When to Use BDF Methods

BDF methods are recommended for:

  - **Large stiff systems** with more than 1000 equations
  - **Very stiff problems** where other implicit methods struggle
  - **Long-time integration** of stiff systems
  - **Parabolic PDEs** after spatial discretization
  - **Reaction-diffusion systems** and chemical kinetics
  - **Circuit simulation** and other engineering applications with large stiff systems

## Solver Selection Guide

### Recommended methods

  - **`QNDF`**: Adaptive order quasi-constant timestep BDF, best general choice for large systems
  - **`FBDF`**: Fixed-leading coefficient BDF, often more efficient than QNDF

## Performance Characteristics

  - **Most efficient for systems with >1000 equations**
  - **Outperform Runge-Kutta methods** on very large stiff systems
  - **Memory efficient** due to multistep structure
  - **Excel at very low tolerances** (1e-9 and below)
  - **Particularly effective** for problems arising from PDE discretizations

## Comparison with Other Methods

Choose BDF methods over:

  - **Rosenbrock methods**: When system size > 1000 equations
  - **SDIRK methods**: For very large stiff systems where RK methods become expensive
  - **Explicit methods**: For any stiff problem

Choose other methods over BDF when:

  - **System size < 100**: Rosenbrock or SDIRK methods often more efficient
  - **Problems with large complex eigenvalues**: Rosenbrock and L-stable SDIRK methods are more stable due to BDF methods only being alpha-stable
  - **Moderate stiffness**: SDIRK methods may be more robust
  - **Non-stiff problems**: Use explicit methods like Tsit5

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqBDF", "QNDF")
```

## Full list of solvers

```@docs
ABDF2
QNDF
QNDF1
QNDF2
QBDF
QBDF1
QBDF2
MEBDF2
FBDF
```
