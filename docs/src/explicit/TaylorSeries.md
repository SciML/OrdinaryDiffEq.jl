```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqTaylorSeries

Taylor series methods for ordinary differential equations using automatic differentiation. These methods achieve very high-order accuracy by computing Taylor expansions of the solution using automatic differentiation techniques through TaylorDiff.jl.

!!! warn "Development Status"
    
    These methods are still in development and may not be fully optimized or reliable for production use.

## Key Properties

Taylor series methods provide:

  - **Very high-order accuracy** with arbitrary order capability
  - **Automatic differentiation** for derivative computation
  - **Step size control** through Taylor series truncation
  - **Natural error estimation** from higher-order terms
  - **Excellent accuracy** for smooth problems
  - **Single-step methods** without requiring history

## When to Use Taylor Series Methods

These methods are recommended for:

  - **Ultra-high precision problems** where maximum accuracy is needed
  - **Smooth problems** with well-behaved derivatives
  - **Scientific computing** requiring very low error tolerances
  - **Problems with expensive function evaluations** where high-order methods reduce total steps

## Mathematical Background

Taylor series methods compute the solution using Taylor expansions:
`u(t + h) = u(t) + h*u'(t) + h²/2!*u''(t) + h³/3!*u'''(t) + ...`

The derivatives are computed automatically using automatic differentiation, allowing arbitrary-order methods without manual derivative computation.

## Solver Selection Guide

### Available Methods

  - **`ExplicitTaylor2`**: Second-order Taylor series method for moderate accuracy
  - **`ExplicitTaylor`**: Arbitrary-order Taylor series method (specify order with `order = Val{p}()`)

### Usage considerations

  - **Smooth problems only**: Methods assume the function has many continuous derivatives
  - **Computational cost**: Higher orders require more automatic differentiation computations
  - **Memory requirements**: Higher orders store more derivative information

## Performance Guidelines

### When Taylor series methods excel

  - **Very smooth problems** where high-order derivatives exist and are well-behaved
  - **High precision requirements** beyond standard double precision
  - **Long-time integration** where accumulated error matters
  - **Problems where function evaluations dominate** computational cost

### Problem characteristics

  - **Polynomial and analytic functions** work extremely well
  - **Smooth ODEs** from physics simulations
  - **Problems requiring** very low tolerances (< 1e-12)

## Limitations and Considerations

### Method limitations

  - **Requires smooth functions** - non-smooth problems may cause issues
  - **Memory overhead** for storing multiple derivatives
  - **Limited to problems** where high-order derivatives are meaningful
  - **Automatic differentiation compatibility** - requires functions compatible with TaylorDiff.jl and Symbolics.jl tracing
  - **Long compile times** due to automatic differentiation and symbolic processing overhead

### When to consider alternatives

  - **Non-smooth problems**: Use adaptive Runge-Kutta methods instead
  - **Stiff problems**: Taylor methods are explicit and may be inefficient
  - **Large systems**: Automatic differentiation cost may become prohibitive
  - **Standard accuracy needs**: Lower-order methods are often sufficient

## Alternative Approaches

Consider these alternatives:

  - **High-order Runge-Kutta** methods (Feagin, Verner) for good accuracy with less overhead
  - **Extrapolation methods** for high accuracy with standard function evaluations
  - **Adaptive methods** for problems with varying smoothness
  - **Implicit methods** for stiff problems requiring high accuracy

## Installation and Usage

Taylor series methods require explicit installation of the specialized library:

```julia
using Pkg
Pkg.add("OrdinaryDiffEqTaylorSeries")
```

Then use the methods with:

```julia
using OrdinaryDiffEqTaylorSeries

# Example: Second-order Taylor method
function f(u, p, t)
    σ, ρ, β = p
    du1 = σ * (u[2] - u[1])
    du2 = u[1] * (ρ - u[3]) - u[2]
    du3 = u[1] * u[2] - β * u[3]
    [du1, du2, du3]
end

u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 10.0)
p = [10.0, 28.0, 8/3]
prob = ODEProblem(f, u0, tspan, p)

# Second-order Taylor method
sol = solve(prob, ExplicitTaylor2())

# Arbitrary-order Taylor method (e.g., 8th order)
sol = solve(prob, ExplicitTaylor(order = Val{8}()))
```

## Full list of solvers

```@docs
ExplicitTaylor2
ExplicitTaylor
```
