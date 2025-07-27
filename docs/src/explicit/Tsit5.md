```@meta
CollapsedDocStrings = true
```

# [OrdinaryDiffEqTsit5](@id OrdinaryDiffEqTsit5)

The Tsitouras 5/4 method is the **recommended default solver** for most non-stiff differential equation problems. This method provides an excellent balance of efficiency, reliability, and accuracy.

## Key Properties

Tsit5 offers:

  - **Fifth-order accuracy** with embedded fourth-order error estimation
  - **Excellent efficiency** at default tolerances (1e-6 to 1e-3)
  - **FSAL (First Same As Last)** property for computational efficiency
  - **High-quality interpolation** for dense output
  - **Robust performance** across a wide range of problem types
  - **Optimized coefficients** for minimal error in practical applications

## When to Use Tsit5

Tsit5 is recommended for:

  - **Most non-stiff problems** as the first choice solver
  - **Default and higher tolerances** (1e-3 to 1e-6)
  - **General-purpose integration** when problem characteristics are unknown
  - **Educational and research applications** as a reliable baseline
  - **Real-time applications** requiring predictable performance
  - **Problems where simplicity and reliability** are preferred over maximum efficiency

## Solver Selection Guide

### Primary recommendation

  - **`Tsit5`**: Use as the default choice for non-stiff problems at standard tolerances

### Automatic switching

  - **`AutoTsit5`**: Automatically switches to a stiff solver when stiffness is detected, making it robust for problems of unknown character

## When to Consider Alternatives

Consider other solvers when:

  - **Higher accuracy needed**: Use Verner methods (Vern6, Vern7, Vern8, Vern9) for tolerances below 1e-6
  - **Higher tolerances**: Use BS3 or OwrenZen3 for tolerances above 1e-3
  - **Robust error control needed**: Use BS5 when Tsit5 struggles with error estimation
  - **Equation is stiff**: Use implicit methods (SDIRK, BDF) or semi-implicit methods (Rosenbrock) for stiff problems
  - **Special properties required**: Use specialized methods (SSP, symplectic, etc.) for specific problem types

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqTsit5", "Tsit5")
```

## Full list of solvers

```@docs
Tsit5
AutoTsit5
```
