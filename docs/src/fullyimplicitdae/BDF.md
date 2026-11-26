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
    out[1] = -0.04u[1] + 1.0e4 * u[2] * u[3] - du[1]
    out[2] = +0.04u[1] - 3.0e7 * u[2]^2 - 1.0e4 * u[2] * u[3] - du[2]
    out[3] = u[1] + u[2] + u[3] - 1.0
    return
end
u₀ = [1.0, 0, 0]
du₀ = [-0.04, 0.04, 0.0]
tspan = (0.0, 100000.0)
differential_vars = [true, true, false]
prob = DAEProblem(f2, du₀, u₀, tspan; differential_vars)
sol = solve(prob, DFBDF())
```

## Sparse Jacobians

`DFBDF` and the other BDF DAE methods accept an `ADTypes.AutoSparse` `autodiff`. For in-place
DAEs with AD-computed separated Jacobians in concrete CSC storage, each partial Jacobian has its
own matrix with the same stored pattern. Both differentiation preparations color that pattern at
integrator construction and reuse their colorings.

  - **The pattern must hold every entry seen during the solve.** An automatically detected
    pattern (`TracerSparsityDetector`, `TracerLocalSparsityDetector`) fixes the storage when
    the integrator is built, from the combined stage Jacobian at its initial point. Every fresh
    Jacobian evaluation re-runs the detector on both partial Jacobians at their real
    `(t, u, du)` and throws a `DimensionMismatch` if it finds an entry the storage does not
    hold, rather than returning a wrong Jacobian. If the detector cannot trace the residual at the
    refresh point (for example a comparison on a global tracer), user-supplied storage skips
    validation for that refresh. Automatically derived storage is instead differentiated column by
    column without compression at that point, with each column perturbed independently (for a
    finite-difference backend the step of a column is scaled by that coordinate alone and uses the
    backend's own stencil and step settings). If every
    nonzero lies inside the stored pattern, those validated values are the partial Jacobian that is
    used; otherwise a `DimensionMismatch` names the first missing entry. For a reverse-mode
    backend (for example `AutoEnzyme(mode = Enzyme.Reverse)`) the same fallback streams one output row
    per reverse sweep through reused buffers instead, so its cost grows with the number of rows.
    Storage that already holds
    every entry is never validated. Storage is never grown during the solve. For a
    state-dependent pattern, supply a `KnownJacobianSparsityDetector` holding the union of all
    its patterns (or a `jac_prototype` holding it), or use a dense `autodiff`.
  - **A known pattern is trusted.** With a `KnownJacobianSparsityDetector` the pattern is not
    re-detected during the solve, so it must hold every entry the residual can ever have.
  - **A supplied `jac_prototype` or `sparsity` must cover every entry.** When a supplied
    `jac_prototype` or `sparsity` defines storage, it must cover every entry. A tracer checks it
    at each fresh Jacobian evaluation it can trace. A known pattern is checked against it when the
    integrator is built; either throws a `DimensionMismatch` if the detector finds an entry the
    storage does not hold. Stored zeros of a known pattern count as stored entries.
  - A global tracer (`TracerSparsityDetector`) requires a traceable residual when the integrator
    builds the combined stage Jacobian, including when a `jac_prototype` supplies the storage. Use
    `TracerLocalSparsityDetector` or `KnownJacobianSparsityDetector` if the residual branches on `u`
    or `du`. Tracer validation and any uncompressed fallback evaluations are excluded from
    `stats.nf`. Their cost depends on the detector, backend, and problem size.
  - The number of colors, and so of residual evaluations per Jacobian, depends on the union
    pattern, not on each partial Jacobian alone.

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
