```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqStabilizedRK

Stabilized Runge-Kutta methods are explicit schemes designed to handle moderately stiff problems by extending the stability region through careful tableau construction. These methods use an upper bound on the spectral radius of the Jacobian to achieve much larger stable timesteps than conventional explicit methods. **These methods are good for large real eigenvalue problems, but not for problems where the complex eigenvalues are large.**

## Key Properties

Stabilized RK methods provide:

  - **Extended stability regions** for moderately stiff problems
  - **Explicit formulation** avoiding nonlinear solvers
  - **Large stable timestep sizes** compared to standard explicit methods
  - **Automatic spectral radius estimation** or user-supplied bounds
  - **Efficient for parabolic PDEs** with moderate stiffness
  - **Good performance** on problems with well-separated timescales

## When to Use Stabilized RK Methods

These methods are recommended for:

  - **Moderately stiff problems** where implicit methods are overkill
  - **Parabolic PDEs** with diffusion-dominated behavior
  - **Problems with large spatial grids** where implicit methods become expensive
  - **Systems with well-separated timescales** but not extreme stiffness
  - **Cases where explicit is preferred** but standard methods are unstable
  - **Large-scale problems** where linear algebra cost of implicit methods is prohibitive

## Mathematical Background

Stabilized methods achieve extended stability by constructing tableaus with enlarged stability regions, often using Chebyshev polynomials or orthogonal polynomial techniques. The stable timestep is determined by the spectral radius bound rather than the CFL condition. **Important**: These methods extend stability primarily along the negative real axis, making them effective for large real eigenvalues but ineffective when complex eigenvalues dominate the stiffness.

## Spectral Radius Estimation

Users can supply an upper bound on the spectral radius using:

```julia
eigen_est = (integrator) -> integrator.eigen_est = upper_bound
```

If not provided, the methods include automatic estimation procedures.

## Solver Selection Guide

### Recommended stabilized methods

  - **`ROCK2`**: Second-order ROW-type stabilized method with extended stability
  - **`ROCK4`**: Fourth-order stabilized method for higher accuracy requirements

## Performance Guidelines

### When stabilized methods excel

  - **Large real eigenvalue problems** where stiffness comes from real eigenvalues
  - **Moderate stiffness ratio** (10³ to 10⁶) dominated by real eigenvalues
  - **Large spatial discretizations** where implicit solver cost is high
  - **Very large systems** where stabilized RK methods are more efficient than BDF methods due to no linear algebra requirements
  - **Parabolic PDEs** with diffusion-dominated (real eigenvalue) stiffness
  - **Problems where spectral radius** can be estimated reliably

### When to use alternatives

  - **Complex eigenvalue dominated problems**: Use implicit methods (BDF, SDIRK, Rosenbrock)
  - **Non-stiff problems**: Use standard explicit methods (Tsit5, Verner)

## Usage Considerations

  - **Spectral radius estimation** is crucial for performance
  - **Method efficiency** depends on stiffness ratio
  - **Test against implicit methods** for highly stiff problems
  - **Consider adaptive spectral radius** estimation for varying stiffness

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqStabilizedRK", "ROCK4")
```

## Full list of solvers

```@docs
ROCK2 
ROCK4 
RKC
SERK2
ESERK4
ESERK5
```
