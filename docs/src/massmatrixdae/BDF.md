```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqBDF

BDF (Backward Differentiation Formula) methods for mass matrix differential-algebraic equations (DAEs) and stiff ODEs with singular mass matrices. These methods provide robust, high-order integration for systems with algebraic constraints and mixed differential-algebraic structure.

## Key Properties

Mass matrix BDF methods provide:

  - **DAE capability** for index-1 differential-algebraic equations
  - **Mass matrix support** for singular and non-diagonal mass matrices
  - **High-order accuracy** up to 5th order with good stability
  - **L-stable behavior** for stiff problems with excellent damping
  - **Automatic differentiation** for efficient Jacobian computation
  - **Variable order and stepsize** adaptation for efficiency

## When to Use Mass Matrix BDF Methods

These methods are recommended for:

  - **Differential-algebraic equations (DAEs)** with index-1 structure
  - **Constrained mechanical systems** with holonomic constraints
  - **Electrical circuit simulation** with algebraic loop equations
  - **Chemical reaction networks** with conservation constraints
  - **Multibody dynamics** with kinematic constraints
  - **Semi-explicit DAEs** arising from spatial discretizations

## Mathematical Background

Mass matrix DAEs have the form:
`M du/dt = f(u,t)`

where M is a potentially singular mass matrix. When M is singular, some equations become algebraic constraints rather than differential equations, leading to a DAE system.

## Problem Formulation

Use `ODEFunction` with a `mass_matrix`:

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
sol = solve(prob_mm, FBDF(), reltol = 1e-8, abstol = 1e-8)
```

## Solver Selection Guide

### Recommended Methods

  - **`FBDF`**: **Recommended** - Fixed leading coefficient BDF with excellent stability
  - **`QNDF`**: Quasi-constant stepsize Nordsieck BDF with good efficiency
  - **`QBDF`**: Alternative quasi-constant stepsize BDF formulation

### Specific order methods

  - **`QNDF1`**: First-order method for simple problems
  - **`QNDF2`**: Second-order method balancing accuracy and stability
  - **`QBDF1`**, **`QBDF2`**: Alternative second-order formulations
  - **`ABDF2`**: Adams-type BDF for specific applications
  - **`MEBDF2`**: Modified extended BDF for enhanced stability

## Performance Guidelines

### When mass matrix BDF methods excel

  - **Index-1 DAE systems** with well-separated differential and algebraic variables
  - **Large stiff systems** with algebraic constraints
  - **Problems with conservation laws** naturally expressed as constraints
  - **Multiphysics simulations** combining differential and algebraic equations
  - **Systems where constraints are essential** to the physics

### Mass matrix considerations

  - **Singular mass matrices** require consistent initial conditions
  - **Index determination** affects solver performance and stability
  - **Constraint violations** may accumulate and require projection
  - **Well-conditioned problems** generally perform better

## Important Considerations

### Initial conditions

  - **Must be consistent** with algebraic constraints
  - **Use initialization procedures** if constraints are not satisfied initially
  - **Index-1 assumption** requires that constraints uniquely determine algebraic variables

### Numerical challenges

  - **Constraint drift** may occur over long integrations
  - **Index higher than 1** not directly supported
  - **Ill-conditioned mass matrices** can cause numerical difficulties
  - **Discontinuities** in constraints require special handling

## Alternative Approaches

Consider these alternatives:

  - **Implicit Runge-Kutta methods** for higher accuracy requirements
  - **Rosenbrock methods** for moderately stiff DAEs
  - **Projection methods** for constraint preservation
  - **Index reduction techniques** for higher-index DAEs

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqBDF", "FBDF")
```

## Full list of solvers

```@docs; canonical=false
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
