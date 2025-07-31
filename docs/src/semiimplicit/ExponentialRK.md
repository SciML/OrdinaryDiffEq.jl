```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqExponentialRK

Exponential Runge-Kutta methods are specialized integrators for semi-linear differential equations of the form `du/dt = Au + f(u,t)`, where `A` is a linear operator (often representing diffusion or dispersion) and `f` represents nonlinear terms. These methods are particularly effective for stiff linear parts combined with non-stiff nonlinear terms. **Important**: The nonlinear term `f(u,t)` must be non-stiff for these methods to be effective.

## Key Properties

Exponential RK methods provide:

  - **Exact treatment of linear parts** using matrix exponential functions
  - **High-order accuracy** for both linear and nonlinear components
  - **Excellent stability properties** for problems with stiff linear operators
  - **Efficient handling of semi-linear PDEs** after spatial discretization
  - **Reduced timestep restrictions** compared to traditional explicit methods
  - **Preservation of qualitative behavior** for many physical systems

## When to Use Exponential RK Methods

These methods are recommended for:

  - **Semi-linear PDEs** with stiff diffusion/dispersion and moderate **non-stiff** nonlinearity
  - **Reaction-diffusion systems** with fast diffusion and slower **non-stiff** reactions
  - **Nonlinear Schrödinger equations** and other dispersive wave equations with **non-stiff** nonlinear terms
  - **Pattern formation problems** (Turing patterns, phase field models) where nonlinearity is **non-stiff**
  - **Quantum dynamics** with linear Hamiltonian and **non-stiff** nonlinear interactions
  - **Problems with strong linear damping** or oscillatory linear parts combined with **non-stiff** nonlinear terms
  - **Spatially discretized PDEs** where the linear part dominates stiffness but the nonlinear part remains **non-stiff**

## Mathematical Background

For problems `du/dt = Au + f(u,t)`, exponential methods compute the exact solution of the linear part `Au` using `exp(A*dt)` and treat the nonlinear part `f(u,t)` with Runge-Kutta-like stages. This approach is particularly effective when `A` represents well-understood physics (diffusion, dispersion, linear oscillations).

## Solver Selection Guide

### Basic exponential time differencing (ETD)

  - **`LawsonEuler`**: First-order exponential Euler method
  - **`NorsettEuler`** / **`ETD1`**: Alternative first-order scheme
  - **`ETDRK2`**: Second-order exponential RK
  - **`ETDRK3`**: Third-order exponential RK
  - **`ETDRK4`**: Fourth-order exponential RK, popular choice
  - **`ETD2`**: Second-order exponential time differencing (in development)

### High-order specialized methods

  - **`HochOst4`**: Fourth-order exponential RK with enhanced stability
  - **`Exp4`**: Fourth-order EPIRK scheme

### Adaptive exponential Rosenbrock

  - **`Exprb32`**: Third-order adaptive method with error control
  - **`Exprb43`**: Fourth-order adaptive method

### EPIRK (Exponential Propagation Iterative RK) methods

  - **`EPIRK4s3A`**: Fourth-order with stiff order 4
  - **`EPIRK4s3B`**: Alternative fourth-order variant
  - **`EPIRK5s3`**: Fifth-order method (note: marked as broken)
  - **`EXPRB53s3`**: Fifth-order with stiff order 5
  - **`EPIRK5P1`**, **`EPIRK5P2`**: Fifth-order variants

## Performance Recommendations

  - **For most semi-linear problems**: `ETDRK4`
  - **For adaptive stepsize**: `Exprb43`
  - **For high stiffness in linear part**: `EPIRK4s3A` or `EPIRK4s3B`
  - **For maximum accuracy**: `EXPRB53s3`

## Implementation Requirements

These methods require:

  - **Computation of matrix exponentials** `exp(A*dt)` and related functions
  - **Krylov subspace methods** for large systems (automatic in most cases)
  - **Proper problem formulation** with identified linear and nonlinear parts

## Installation

To be able to access the solvers in `OrdinaryDiffEqLinear`, you must first install them use the Julia package manager:

```julia
using Pkg
Pkg.add("OrdinaryDiffEqExponentialRK")
```

This will only install the solvers listed at the bottom of this page.
If you want to explore other solvers for your problem,
you will need to install some of the other libraries listed in the navigation bar on the left.

## Example usage

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqExponentialRK", "EPIRK5s3")
```

## Full list of solvers

```@docs; canonical=false
LawsonEuler
NorsettEuler
ETD2
ETDRK2
ETDRK3
ETDRK4
HochOst4
```

### Adaptive Exponential Rosenbrock Methods

```@docs
Exprb32
Exprb43
```

### Exponential Propagation Iterative Runge-Kutta Methods (EPIRK)

```@docs
Exp4
EPIRK4s3A
EPIRK4s3B
EPIRK5s3
EXPRB53s3
EPIRK5P1
EPIRK5P2
```
