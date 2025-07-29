```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqAdamsBashforthMoulton

Adams-Bashforth and Adams-Moulton multistep methods for non-stiff differential equations. **Note that Runge-Kutta methods generally come out as more efficient in benchmarks, except when the ODE function `f` is expensive to evaluate or the problem is very smooth.** These methods can achieve high accuracy with fewer function evaluations per step than Runge-Kutta methods in those specific cases.

## Key Properties

Adams-Bashforth-Moulton methods provide:

  - **Reduced function evaluations** compared to Runge-Kutta methods
  - **High efficiency** for expensive-to-evaluate functions
  - **Multistep structure** using information from previous timesteps
  - **Variable step and order** capabilities for adaptive integration
  - **Predictor-corrector variants** for enhanced accuracy and stability
  - **Good stability properties** for non-stiff problems

## When to Use Adams-Bashforth-Moulton Methods

These methods are recommended for:

  - **Expensive function evaluations** where minimizing calls to `f` is critical
  - **Non-stiff smooth problems** with regular solution behavior
  - **Long-time integration** where efficiency over many steps matters
  - **Problems with expensive Jacobian computations** that cannot use implicit methods efficiently
  - **Scientific computing applications** with computationally intensive right-hand sides
  - **Systems where startup cost** of multistep methods is amortized over long integration

## Method Types

### Explicit Adams-Bashforth (AB)

Pure explicit multistep methods using only past information:

  - **Lower computational cost** per step
  - **Less stability** than predictor-corrector variants
  - **Good for mildly stiff** problems

### Predictor-Corrector Adams-Bashforth-Moulton (ABM)

Implicit corrector step for enhanced accuracy:

  - **Better accuracy** than pure explicit methods
  - **Improved stability** properties
  - **Slightly higher cost** but often worth it

## Solver Selection Guide

### Primary recommendation

  - **`VCABM`**: **Main recommendation** - adaptive order variable-step Adams-Bashforth-Moulton, best overall choice for Adams methods

### Variable-step predictor-corrector methods

  - **`VCABM3`**: Third-order variable-step Adams-Bashforth-Moulton
  - **`VCABM4`**: Fourth-order variable-step Adams-Bashforth-Moulton
  - **`VCABM5`**: Fifth-order variable-step Adams-Bashforth-Moulton

### Variable-step Adams-Bashforth methods

  - **`VCAB3`**: Third-order variable-step Adams-Bashforth
  - **`VCAB4`**: Fourth-order variable-step Adams-Bashforth
  - **`VCAB5`**: Fifth-order variable-step Adams-Bashforth

### Fixed-step predictor-corrector methods

  - **`ABM32`**: Third-order Adams-Bashforth-Moulton
  - **`ABM43`**: Fourth-order Adams-Bashforth-Moulton
  - **`ABM54`**: Fifth-order Adams-Bashforth-Moulton

### Fixed-step explicit methods

  - **`AB3`**: Third-order Adams-Bashforth
  - **`AB4`**: Fourth-order Adams-Bashforth
  - **`AB5`**: Fifth-order Adams-Bashforth

## Performance Considerations

  - **Most efficient** when function evaluation dominates computational cost
  - **Startup phase** requires initial steps from single-step method
  - **Memory efficient** compared to high-order Runge-Kutta methods
  - **Best for smooth problems** - avoid for problems with discontinuities

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqAdamsBashforthMoulton", "VCABM")
```

## Full list of solvers

### Explicit Multistep Methods

```@docs
AB3
AB4
AB5
```

### Predictor-Corrector Methods

```@docs
ABM32
ABM43
ABM54
VCAB3
VCAB4
VCAB5
VCABM3
VCABM4
VCABM5
VCABM
```
