```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqHighOrderRK

High-order explicit Runge-Kutta methods for non-stiff differential equations requiring high accuracy. These methods provide alternatives to the Verner methods, though **[`OrdinaryDiffEqVerner`](@ref OrdinaryDiffEqVerner) methods generally perform better at low tolerances** and should be preferred in most cases.

## Key Properties

High-order RK methods provide:

  - **High-order accuracy** (7th and 8th order) for precise integration
  - **Specialized coefficients** for specific problem types
  - **Dense output capabilities** for some methods
  - **Alternative approaches** to the more commonly used Verner methods

## When to Use High-Order RK Methods

These methods are recommended when:

  - **Verner methods are not suitable** for specific problem characteristics
  - **Specialized properties** are needed (e.g., phase-fitted methods for oscillatory problems)
  - **Research or comparison purposes** require different high-order method families
  - **Specific applications** benefit from particular method properties

## Solver Selection Guide

### General high-order integration

  - **Use [`Vern7`](@ref OrdinaryDiffEqVerner) or [`Vern8`](@ref OrdinaryDiffEqVerner) instead** - they are generally more efficient

### Specialized cases where these methods may be preferred

  - **`TanYam7`**: Seventh-order Tanaka-Yamashita method
  - **`TsitPap8`**: Eighth-order Tsitouras-Papakostas method
  - **`DP8`**: Eighth-order Dormand-Prince method (Hairer's 8/5/3 implementation)
  - **`PFRK87`**: Phase-fitted eighth-order method for oscillatory problems

## Performance Notes

  - **Verner methods are generally more efficient** for most high-accuracy applications
  - **These methods are provided** for specialized use cases and research purposes
  - **Consider problem-specific properties** when choosing between different high-order families

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqHighOrderRK", "DP8")
```

## Full list of solvers

```@docs
TanYam7
TsitPap8
DP8
PFRK87
```
