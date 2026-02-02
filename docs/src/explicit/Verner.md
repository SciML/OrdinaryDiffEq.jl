```@meta
CollapsedDocStrings = true
```

# [OrdinaryDiffEqVerner](@id OrdinaryDiffEqVerner)

Verner methods are high-order explicit Runge-Kutta methods designed for high-accuracy integration of non-stiff differential equations. These are the preferred solvers when very low tolerances are required.

## Key Properties

Verner methods provide:

  - **High-order accuracy** (6th through 9th order) for precise integration
  - **Excellent efficiency** at low tolerances (1e-8 to 1e-15)
  - **Robust error estimation** with embedded error control
  - **Dense output capability** with high-quality interpolation

## When to Use Verner Methods

Verner methods are recommended for:

  - **High-accuracy requirements** with tolerances between 1e-8 and 1e-15
  - **Smooth non-stiff problems** where high precision is critical
  - **Long-time integration** where error accumulation must be minimized
  - **Problems requiring dense output** with high interpolation accuracy
  - **Orbit computation, molecular dynamics,** and other precision-critical applications

## Solver Selection Guide

### Medium-low tolerance (1e-6 to 1e-8)

  - **`Vern6`**: Sixth-order method, good balance of efficiency and accuracy
  - **`AutoVern6`**: Automatic switching version for mixed stiffness

### Low tolerance (1e-8 to 1e-12) with Float64

  - **`Vern7`**: Seventh-order method, excellent for most high-precision needs
  - **`Vern8`**: Eighth-order method, best efficiency at very low tolerances
  - **`AutoVern7`**, **`AutoVern8`**: Automatic switching versions

### Very low tolerance (<1e-12)

  - **`Vern9`**: Ninth-order method for extreme precision requirements
  - **Recommended with `BigFloat`** for tolerances below 1e-15
  - **`AutoVern9`**: Automatic switching version for mixed problems

## Performance Notes

  - **Vern6**: Most efficient for tolerances around 1e-6 to 1e-8
  - **Vern7**: Sweet spot for tolerances around 1e-8 to 1e-10
  - **Vern8**: Best for tolerances around 1e-10 to 1e-12
  - **Vern9**: For tolerances below 1e-12, especially with arbitrary precision

The `Auto*` variants automatically switch to stiff solvers when stiffness is detected, making them robust for problems of unknown character.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqVerner", "Vern6")
```

## Full list of solvers

```@docs
Vern6
Vern7
Vern8
Vern9
```

```@docs
AutoVern6
AutoVern7
AutoVern8
AutoVern9
```
