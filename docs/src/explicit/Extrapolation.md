```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqExtrapolation

Explicit extrapolation methods that achieve high accuracy through Richardson extrapolation of basic integration schemes. These methods provide adaptive order capabilities and natural parallelism, though they are generally outclassed by modern Runge-Kutta methods for most non-stiff problems.

## Key Properties

Extrapolation methods provide:

  - **Adaptive order capability** allowing arbitrarily high orders
  - **Natural parallelism** across different substep sequences
  - **High accuracy potential** for very smooth problems
  - **Richardson extrapolation** to eliminate lower-order error terms
  - **Automatic stepsize and order control**
  - **Theoretical appeal** but often practical limitations

## When to Use Extrapolation Methods

These methods are recommended for:

  - **Very smooth problems** where high-order accuracy is beneficial
  - **Extremely low tolerance requirements** where adaptive order helps
  - **Parallel computing environments** that can exploit the natural parallelism
  - **Research applications** exploring adaptive order techniques
  - **Problems where other high-order methods struggle** with accuracy

## Important Limitations

  - **Generally outclassed** by modern explicit RK methods (Tsit5, Verner methods)
  - **Higher computational overhead** compared to optimized RK methods
  - **Best suited for very smooth functions** - poor performance on non-smooth problems
  - **Parallel efficiency gains** often don't compensate for increased work

## Mathematical Background

Extrapolation methods use sequences of basic integrators (like Euler or midpoint) with different stepsizes, then apply Richardson extrapolation to achieve higher-order accuracy. The adaptive order capability comes from using longer extrapolation sequences.

## Solver Selection Guide

### Explicit extrapolation methods

  - **`AitkenNeville`**: Euler extrapolation using Aitken-Neville algorithm
  - **`ExtrapolationMidpointDeuflhard`**: Midpoint extrapolation with barycentric coordinates
  - **`ExtrapolationMidpointHairerWanner`**: Midpoint extrapolation following ODEX algorithm

### When to consider these methods

  - **Very low tolerances** (< 1e-12) where adaptive order might help
  - **Extremely smooth problems** with analytic solutions
  - **Parallel computing** scenarios with many available cores
  - **Comparison studies** with other high-order methods

### Better alternatives for most problems

  - **For high accuracy**: Use Verner methods (Vern7, Vern8, Vern9)
  - **For general problems**: Use Tsit5 or appropriate RK method
  - **For stiff problems**: Consider [implicit extrapolation methods](@ref StiffExtrapolation)

## Performance Notes

  - **Consider stiff extrapolation** methods which can perform very well for sufficiently stiff problems
  - **Test against Verner methods** before choosing extrapolation for high accuracy
  - **Parallelism benefits** are problem and hardware dependent
  - **Most effective** on very smooth, well-behaved problems

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqExtrapolation", "ExtrapolationMidpointDeuflhard")
```

## Full list of solvers

```@docs
AitkenNeville
ExtrapolationMidpointDeuflhard
ExtrapolationMidpointHairerWanner
```
