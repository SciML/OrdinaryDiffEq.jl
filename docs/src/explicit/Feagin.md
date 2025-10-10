```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqFeagin

Ultra-high-order explicit Runge-Kutta methods for non-stiff problems at extremely low tolerances (< 1e-30). These methods are designed for applications requiring extreme precision, typically used with higher-precision number types like `BigFloat`.

## Key Properties

Feagin methods provide:

  - **Ultra-high-order accuracy** (10th, 12th, and 14th order)
  - **Extreme precision capabilities** for very low tolerance requirements
  - **Compatibility with arbitrary precision** arithmetic (`BigFloat`, `Float128`)
  - **Specialized for very demanding applications** requiring maximum accuracy

## When to Use Feagin Methods

These methods are recommended for:

  - **Extremely low tolerance problems** (< 1e-30)
  - **Arbitrary precision arithmetic** applications using `BigFloat` or `Float128`
  - **Ultra-high precision requirements** where standard methods are insufficient
  - **Research applications** requiring maximum possible accuracy
  - **Long-time integration** where error accumulation must be minimized to extreme levels

## Important Limitations

### Theoretical vs Practical Performance

  - **Very good theoretical efficiency** due to high order and optimized coefficients
  - **Poor practical performance** in benchmarks due to bad error estimators and adaptivity issues
  - **Generally recommend `Vern9` instead** as it tends to be more efficient in practice despite lower theoretical order

### Performance Considerations

  - **May be less efficient** than `Vern9` even for very low tolerance problems
  - **Outperformed by extrapolation methods** at extremely low tolerances due to adaptive order
  - **Potential efficiency for >128-bit numbers** but no practical cases found yet where this is actually true
  - **Should always be tested** against `Vern9` and extrapolation methods

## Solver Selection Guide

### Extreme precision (< 1e-30)

  - **`Feagin14`**: 14th-order method for maximum accuracy
  - **`Feagin12`**: 12th-order method, balance of accuracy and efficiency
  - **`Feagin10`**: 10th-order method for moderate extreme precision

### Strongly recommended alternatives

  - **For most very low tolerance problems**: Use `Vern9` first (more efficient in practice despite lower theoretical order)
  - **For extremely low tolerances**: Consider extrapolation methods for adaptive order
  - **For >128-bit precision**: These methods may be more efficient, but no practical cases found yet
  - **Always benchmark**: Compare performance with `Vern9` and extrapolation methods before choosing Feagin methods

## Usage Guidelines

  - **Best with `BigFloat`** or `Float128` number types
  - **Useful in Float128 precision range** but test against other algorithms
  - **Consider problem-specific characteristics** when choosing order level

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqFeagin", "Feagin14")
```

## Full list of solvers

```@docs
Feagin10
Feagin12
Feagin14
```
