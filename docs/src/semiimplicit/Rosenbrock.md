```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqRosenbrock

Rosenbrock methods are semi-implicit Runge-Kutta methods designed for small to medium-sized stiff systems of differential equations. These methods combine the stability properties of implicit methods with computational efficiency by using approximate Jacobians. They work particularly well for strict tolerance requirements.

## Key Properties

Rosenbrock methods provide:

  - **Excellent efficiency for small to medium systems** (<1000 ODEs)
  - **L-stable and A-stable variants** for stiff problems
  - **W-method structure** making them robust to inaccurate Jacobians
  - **Automatic differentiation compatibility** for Jacobian computation
  - **High-order accuracy** with embedded error estimation
  - **Excellent performance for strict tolerance requirements**
  - **Stiffly accurate variants** for enhanced stability

## When to Use Rosenbrock Methods

Rosenbrock methods are recommended for:

  - **Small to medium stiff systems** (10 to 1000 equations)
  - **Problems where Jacobians are available** or can be computed efficiently
  - **Medium tolerance requirements** (1e-8 to 1e-2)
  - **Stiff ODEs arising from reaction kinetics** and chemical systems
  - **Moderately stiff PDEs** after spatial discretization
  - **Problems requiring reliable error control** for stiff systems

## Solver Selection Guide

### Low tolerance (>1e-2)

  - **`Rosenbrock23`**: Second/third-order method, recommended for low tolerance requirements

### Medium tolerance (1e-8 to 1e-2)

  - **`Rodas5P`**: Fifth-order method, most efficient for many problems
  - **`Rodas4P`**: Fourth-order method, more reliable than Rodas5P
  - **`Rodas5Pe`**: Enhanced fifth-order variant with stiffly accurate embedded estimate for better adaptivity on highly stiff equations
  - **`Rodas5Pr`**: Fifth-order variant with residual test for robust error estimation that guarantees accuracy on interpolation

## Performance Guidelines

  - **Rodas5P**: Best overall efficiency at medium tolerances
  - **Rodas4P**: Most reliable when Rodas5P fails
  - **Rosenbrock23**: Fastest at high tolerances (>1e-2)

## When to Choose Alternatives

Consider other methods when:

  - **System size > 1000**: Use BDF methods (QNDF, FBDF)
  - **Matrix-free methods**: If using Krylov solvers for the linear solver (matrix-free methods), SDIRK or BDF methods are preferred

## Advantages of Rosenbrock Methods

  - **Very low tolerances**: Rosenbrock23 performs well even at very low tolerances
  - **Very high stiffness**: Rosenbrock methods are often more stable than SDIRK or BDF methods because other methods can diverge due to bad initial guesses for the Newton method (leading to nonlinear solver divergence)

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqRosenbrock", "Rodas5P")
```

## Full list of solvers

```@docs
Rosenbrock23
Rosenbrock32
ROS3P
Rodas3
Rodas23W
Rodas3P
Rodas4
Rodas42
Rodas4P
Rodas4P2
Rodas5
Rodas5P
Rodas5Pe
Rodas5Pr
RosenbrockW6S4OS
ROS2
ROS2PR
ROS2S
ROS3
ROS3PR
Scholz4_7
ROS34PW1a
ROS34PW1b
ROS34PW2
ROS34PW3
ROS34PRw
ROS3PRL
ROS3PRL2
ROK4a
RosShamp4
Veldd4
Velds4
GRK4T
GRK4A
Ros4LStab
```
