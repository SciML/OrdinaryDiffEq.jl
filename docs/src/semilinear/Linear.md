```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqLinear

Specialized methods for linear and semi-linear differential equations where the system can be written in matrix form `du/dt = A(t,u) * u` or `du/dt = A(t) * u`. These methods exploit the linear structure to provide exact solutions or highly accurate integration.

## Key Properties

Linear ODE methods provide:

  - **Exact solutions** for time-independent linear systems
  - **Geometric integration** preserving Lie group structure
  - **High-order Magnus expansions** for time-dependent linear systems
  - **Lie group methods** for matrix differential equations
  - **Excellent stability** for a wide range of linear systems
  - **Specialized algorithms** for different types of linear operators

## When to Use Linear Methods

These methods are essential for:

  - **Linear systems** `du/dt = A * u` with constant or time-dependent matrices
  - **Matrix differential equations** on Lie groups (rotation matrices, etc.)
  - **Quantum dynamics** with Hamiltonian evolution
  - **Linear oscillators** and harmonic systems
  - **Time-dependent linear systems** with periodic or smooth coefficients
  - **Geometric mechanics** requiring preservation of group structure

## Mathematical Background

For linear systems `du/dt = A(t) * u`, the exact solution is `u(t) = exp(∫A(s)ds) * u₀`. Linear methods approximate this matrix exponential using various mathematical techniques like Magnus expansions, Lie group integrators, and specialized exponential methods.

## Solver Selection Guide

### Time and state-independent (constant A)

  - **`LinearExponential`**: Exact solution for `du/dt = A * u` with constant matrix A

### Time-dependent, state-independent (A(t))

  - **`MagnusMidpoint`**: Second-order Magnus method
  - **`MagnusLeapfrog`**: Second-order Magnus leapfrog scheme
  - **`MagnusGauss4`**: Fourth-order with Gauss quadrature
  - **`MagnusGL4`**: Fourth-order Gauss-Legendre Magnus method
  - **`MagnusGL6`**: Sixth-order Gauss-Legendre Magnus method
  - **`MagnusGL8`**: Eighth-order Gauss-Legendre Magnus method
  - **`MagnusNC6`**: Sixth-order Newton-Cotes Magnus method
  - **`MagnusNC8`**: Eighth-order Newton-Cotes Magnus method

### State-dependent (A(u))

  - **`LieEuler`**: First-order Lie group method
  - **`RKMK2`**: Second-order Runge-Kutta-Munthe-Kaas method
  - **`RKMK4`**: Fourth-order Runge-Kutta-Munthe-Kaas method
  - **`LieRK4`**: Fourth-order Lie Runge-Kutta method
  - **`CG2`**: Second-order Crouch-Grossman method
  - **`CG4a`**: Fourth-order Crouch-Grossman method
  - **`CayleyEuler`**: First-order method using Cayley transformations

### Adaptive methods

  - **`MagnusAdapt4`**: Fourth-order adaptive Magnus method

### Time and state-dependent (A(t,u))

  - **`CG3`**: Third-order Crouch-Grossman method for most general case

## Method Selection Guidelines

  - **For constant linear systems**: `LinearExponential` (exact)
  - **For time-dependent systems**: Magnus methods based on desired order
  - **For matrix Lie groups**: Lie group methods (RKMK, LieRK4, CG)
  - **For high accuracy**: Higher-order Magnus methods (GL6, GL8)
  - **For adaptive integration**: `MagnusAdapt4`

## Special Considerations

These methods require:

  - **Proper problem formulation** with identified linear structure
  - **Matrix operator interface** for operator-based problems
  - **Understanding of Lie group structure** for geometric problems

## Installation

To be able to access the solvers in `OrdinaryDiffEqLinear`, you must first install them use the Julia package manager:

```julia
using Pkg
Pkg.add("OrdinaryDiffEqLinear")
```

This will only install the solvers listed at the bottom of this page.
If you want to explore other solvers for your problem,
you will need to install some of the other libraries listed in the navigation bar on the left.

## Example usage

```julia
using OrdinaryDiffEqLinear, SciMLOperators
function update_func!(A, u, p, t)
    A[1, 1] = 0
    A[2, 1] = sin(u[1])
    A[1, 2] = -1
    A[2, 2] = 0
end
A0 = ones(2, 2)
A = MatrixOperator(A0, update_func! = update_func!)
u0 = ones(2)
tspan = (0.0, 30.0)
prob = ODEProblem(A, u0, tspan)
sol = solve(prob, LieRK4(), dt = 1 / 4)
```

## Full list of solvers

### Time and State-Independent Solvers

```@docs
LinearExponential
```

### Time-Dependent and State-Independent Solvers

```@docs
MagnusMidpoint
MagnusLeapfrog
MagnusGauss4
MagnusNC6
MagnusGL6
MagnusGL8
MagnusNC8
MagnusGL4
```

### State-Dependent Solvers

```@docs
LieEuler
RKMK2
RKMK4
LieRK4
CG2
CG4a
MagnusAdapt4
CayleyEuler
```

### Time and State-Dependent Operators

```@docs
CG3
```
