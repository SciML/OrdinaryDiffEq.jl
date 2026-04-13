```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqRosenbrock

Rosenbrock methods for mass matrix differential-algebraic equations (DAEs) and stiff ODEs with singular mass matrices. These methods provide efficient integration for moderately stiff systems with algebraic constraints, offering excellent performance for small to medium-sized DAE problems.

## Key Properties

Mass matrix Rosenbrock methods provide:

  - **DAE capability** for index-1 differential-algebraic equations
  - **W-method efficiency** using approximate Jacobians for computational savings
  - **Mass matrix support** for singular and non-diagonal mass matrices
  - **Moderate to high order accuracy** (2nd to 6th order available)
  - **Good stability properties** with stiffly accurate behavior
  - **Embedded error estimation** for adaptive timestepping

## When to Use Mass Matrix Rosenbrock Methods

These methods are recommended for:

  - **Index-1 DAE systems** with moderate stiffness
  - **Small to medium constrained systems** (< 1000 equations)
  - **Semi-explicit DAEs** arising from discretized PDEs
  - **Problems requiring good accuracy** with moderate computational cost
  - **DAEs with moderate nonlinearity** where W-methods are efficient
  - **Electrical circuits** and **mechanical systems** with constraints

!!! warn
    
    In order to use OrdinaryDiffEqRosenbrock with DAEs that require a non-trivial
    consistent initialization, a nonlinear solver is required and thus
    `using OrdinaryDiffEqNonlinearSolve` is required or you must pass an `initializealg`
    with a valid `nlsolve` choice.

## Mathematical Background

Mass matrix DAEs have the form:
`M du/dt = f(u,t)`

Rosenbrock methods linearize around the current solution and solve linear systems of the form:
`(M/γh - J) k_i = ...`

where J is the Jacobian of f and γ is a method parameter.

## Solver Selection Guide

### Recommended Methods by Tolerance

  - **High tolerances** (>1e-2): **`Rosenbrock23`** - efficient low-order method
  - **Medium tolerances** (1e-8 to 1e-2): **`Rodas5P`** - most efficient choice, or **`Rodas4P`** for higher reliability
  - **Low tolerances** (<1e-8): **`Rodas5Pe`** or higher-order alternatives

### Method families

  - **`Rodas5P`**: **Recommended** - Most efficient 5th-order method for general use
  - **`Rodas4P`**: More reliable 4th-order alternative
  - **`Rosenbrock23`**: Good for high tolerance problems
  - **`Rodas5`**: Standard 5th-order method without embedded pair optimization

## Performance Guidelines

### When mass matrix Rosenbrock methods excel

  - **Small to medium DAE systems** (< 1000 equations)
  - **Moderately stiff problems** where full BDF methods are overkill
  - **Problems with efficient Jacobian computation** or finite difference approximation
  - **Index-1 DAEs with well-conditioned mass matrices**
  - **Semi-explicit index-1 problems** from spatial discretizations

### System size considerations

  - **Small systems** (< 100): Rosenbrock methods often outperform multistep methods
  - **Medium systems** (100-1000): Good performance with proper linear algebra
  - **Large systems** (> 1000): Consider BDF methods instead

## Important DAE Considerations

### Initial conditions

  - **Must be consistent** with algebraic constraints
  - **Consistent initialization** may require nonlinear solver
  - **Index-1 assumption** for reliable performance

### Mass matrix requirements

  - **Index-1 DAE structure** for optimal performance
  - **Non-singular leading submatrix** for differential variables
  - **Well-conditioned constraint equations**

## Alternative Approaches

Consider these alternatives:

  - **Mass matrix BDF methods** for larger or highly stiff DAE systems
  - **Implicit Runge-Kutta methods** for higher accuracy requirements
  - **Standard Rosenbrock methods** for regular ODEs without constraints
  - **IMEX methods** if natural explicit/implicit splitting exists

## Example Usage

```julia
using LinearAlgebra: Diagonal
function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    nothing
end
M = Diagonal([1.0, 1.0, 0])  # Singular mass matrix
f = ODEFunction(rober, mass_matrix = M)
prob_mm = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 1e5), (0.04, 3e7, 1e4))
sol = solve(prob_mm, Rodas5(), reltol = 1e-8, abstol = 1e-8)
```

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqRosenbrock", "Rodas5P")
```

## Full list of solvers

```@docs; canonical=false
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
