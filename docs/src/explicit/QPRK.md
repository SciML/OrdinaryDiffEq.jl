```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqQPRK

Quadruple-precision parallel Runge-Kutta (QPRK) methods are high-order explicit solvers specifically designed for ultra-high precision computations using quad-precision arithmetic (`Float128`). These methods combine parallel evaluation capabilities with coefficients optimized for extended precision arithmetic. **Note: These methods are still under-benchmarked and need more research.**

## Key Properties

QPRK methods provide:

  - **Ultra-high-order accuracy** (9th order) for maximum precision
  - **Quadruple-precision optimization** specifically designed for `Float128`
  - **Parallel function evaluations** for computational efficiency
  - **Extreme precision capabilities** for very demanding applications
  - **Optimized coefficients** for extended precision arithmetic

## When to Use QPRK Methods

These methods are recommended for:

  - **Ultra-high precision requirements** demanding `Float128` arithmetic
  - **Extremely low tolerances** (< 1e-20) where standard precision fails
  - **Scientific applications** requiring maximum possible accuracy
  - **Parallel computing environments** with quad-precision support
  - **Research applications** exploring limits of numerical precision
  - **Long-time integration** where error accumulation must be minimized to extreme levels

## Important Requirements

### Precision Requirements

  - **Must use `Float128`** or higher precision number types
  - **All problem components** should support extended precision
  - **Tolerances should match** the precision capabilities (< 1e-20)

### Computational Considerations

  - **Slower** than standard precision methods due to extended precision arithmetic
  - **Higher memory usage** due to extended precision
  - **Limited hardware support** for quad-precision operations

## Mathematical Background

QPRK methods use tableaus with coefficients computed in extended precision to maintain accuracy throughout the ultra-high precision computation. The parallel structure allows independent function evaluations to be computed simultaneously.

## Solver Selection Guide

### Available methods

  - **`QPRK98`**: Ninth-order method optimized for quad-precision arithmetic with parallel evaluation

### Usage guidelines

  - **Essential to use `Float128`** for the state vector and parameters
  - **Consider [MultiFloats.jl](https://github.com/dzhang314/MultiFloats.jl)** for higher precision number types
  - **Set very low tolerances** (e.g., 1e-25) to utilize full precision
  - **Test against alternatives** like Feagin methods with `BigFloat`

## Performance Considerations

  - **Slower** than standard precision methods due to extended precision arithmetic
  - **Memory intensive** due to extended precision storage
  - **Hardware dependent** - some architectures lack efficient quad-precision support

## Alternative High-Precision Methods

For ultra-high precision, also consider:

  - **Feagin methods** with `BigFloat` for arbitrary precision
  - **Arbitrary precision extrapolation** methods
  - **Verner methods** with `BigFloat` for slightly lower but efficient precision
  - **Taylor series methods** with automatic differentiation for extreme precision

## Usage Example

```julia
using OrdinaryDiffEqQPRK
# Ensure using Float128 for ultra-high precision
u0 = Float128[1.0, 0.0]
tspan = (Float128(0.0), Float128(10.0))
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, QPRK98(), abstol = 1e-25, reltol = 1e-25)
```

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqQPRK", "QPRK98")
```

## Full list of solvers

```@docs
QPRK98
```
