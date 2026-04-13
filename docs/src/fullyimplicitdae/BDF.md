```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqBDF

BDF (Backward Differentiation Formula) methods for fully implicit differential-algebraic equations (DAEs) in the form F(du/dt, u, t) = 0. These methods provide robust integration for index-1 DAE systems with fully implicit formulations.

!!! warn "Performance Consideration"
    
    DFBDF and family have not been made fully efficient yet, and thus Sundials.jl IDA is recommended for production use.

## Key Properties

Fully implicit DAE BDF methods provide:

  - **General DAE capability** for F(du/dt, u, t) = 0 formulations
  - **Index-1 DAE support** for properly formulated DAE systems
  - **Robust nonlinear solver integration** for implicit equation systems
  - **High-order accuracy** with excellent stability properties
  - **Large stiff system capability** with efficient linear algebra

## When to Use Fully Implicit DAE BDF Methods

These methods are recommended for:

  - **Fully implicit DAE systems** where F(du/dt, u, t) = 0 cannot be easily rearranged
  - **Index-1 DAE problems** that cannot be easily rearranged to semi-explicit form
  - **Multibody dynamics** with complex kinematic constraints
  - **Electrical circuits** with ideal components and algebraic loops
  - **Chemical engineering** with equilibrium and conservation constraints
  - **Large-scale DAE systems** requiring robust implicit integration

## Mathematical Background

Fully implicit DAEs have the general form:
F(du/dt, u, t) = 0

Unlike semi-explicit forms, these cannot be written as du/dt = f(u,t) even after constraint elimination. BDF methods discretize the time derivative using backward differences and solve the resulting nonlinear system at each timestep.

## Problem Formulation

Use `DAEProblem` with implicit function specification:

```julia
function f2(out, du, u, p, t)
    out[1] = -0.04u[1] + 1e4 * u[2] * u[3] - du[1]
    out[2] = +0.04u[1] - 3e7 * u[2]^2 - 1e4 * u[2] * u[3] - du[2]
    out[3] = u[1] + u[2] + u[3] - 1.0
end
u₀ = [1.0, 0, 0]
du₀ = [-0.04, 0.04, 0.0]
tspan = (0.0, 100000.0)
differential_vars = [true, true, false]
prob = DAEProblem(f2, du₀, u₀, tspan, differential_vars = differential_vars)
sol = solve(prob, DFBDF())
```

## Solver Selection Guide

### Recommended DAE Methods

  - **`DFBDF`**: **Recommended** - Variable-order BDF for general DAE systems
  - **`DImplicitEuler`**: For non-smooth problems with discontinuities

### Method characteristics

  - **`DFBDF`**: Most robust and efficient for general smooth DAE problems
  - **`DImplicitEuler`**: Best choice for problems with discontinuities, events, or non-smooth behavior

## Performance Guidelines

### When fully implicit DAE BDF methods excel

  - **Index-1 DAE systems** with complex implicit structure
  - **Complex constraint structures** with multiple algebraic relationships
  - **Large-scale problems** where specialized DAE methods are essential
  - **Multiphysics simulations** with mixed differential-algebraic structure
  - **Problems where semi-explicit formulation is impractical**

### Index considerations

  - **Index-1 formulation required**: Problems should be written in index-1 form
  - **Compare with mass matrix methods**: For some index-1 problems, mass matrix formulation may be more efficient
  - **Higher-index problems**: Should be reduced to index-1 form before using these methods

## Important DAE Requirements

### Initial conditions

  - **Both u₀ and du₀** must be provided and consistent with constraints
  - **differential_vars** specification helps identify algebraic variables
  - **Consistent initialization** is crucial for index-1 DAE problems

### Function specification

  - **Residual form**: F(du/dt, u, t) = 0 with F returning zero for satisfied equations
  - **Proper scaling**: Ensure equations are well-conditioned numerically
  - **Jacobian availability**: Analytical Jacobians improve performance when available

## Alternative Approaches

Consider these alternatives:

  - **Mass matrix DAE methods** for index-1 problems with M du/dt = f(u,t) structure
  - **Index reduction techniques** using [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) to convert problems to index-1 form if needed
  - **Constraint stabilization** methods for drift control
  - **Projection methods** for manifold preservation

For more details on DAE formulations and alternative approaches, see [this blog post on Neural DAEs](https://www.stochasticlifestyle.com/machine-learning-with-hard-constraints-neural-differential-algebraic-equations-daes-as-a-general-formalism/).

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqBDF", "DFBDF")
```

## Full list of solvers

### DAE

```@docs
DImplicitEuler
DABDF2
DFBDF
```
