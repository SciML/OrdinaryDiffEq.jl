```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqNordsieck

Nordsieck form multistep methods represent an alternative approach to traditional multistep algorithms. Instead of storing past solution values, these methods maintain a vector of scaled derivatives (similar to Taylor series coefficients) to advance the solution. This representation was pioneered in classic codes like LSODE, VODE, and CVODE.

!!! warning "Research and Development"
    
    These methods are currently in research and development and not intended for general use.

## Key Properties

Nordsieck methods provide:

  - **Derivative-based representation** instead of solution history
  - **Improved restartability** after discontinuities using derivative information
  - **Variable order and stepsize** capabilities
  - **Alternative to history-based** multistep methods
  - **Research and experimental** implementations

## When to Use Nordsieck Methods

These methods are recommended for:

  - **Research applications** exploring alternative multistep representations
  - **Problems with discontinuities** where restartability is important
  - **Experimental comparisons** with traditional multistep methods
  - **Development of discontinuity-aware** algorithms

## Important Limitations

### Experimental Status

  - **Considered experimental** and inferior to modern BDF implementations
  - **Generally recommend FBDF instead** for production use
  - **Maintained for research purposes** and future development
  - **Numerical instabilities** can arise from higher derivative representations

### Performance Considerations

  - **Less robust** than fixed-leading coefficient BDF methods
  - **Higher computational overhead** for derivative maintenance
  - **Potential stability issues** with derivative representations

## Mathematical Background

The Nordsieck form represents the solution using scaled derivatives:
`y_n = [y, h*y', h²*y''/2!, h³*y'''/3!, ...]`

This representation allows reconstruction of the solution and its derivatives, enabling restarts after discontinuities without losing accuracy.

## Solver Selection Guide

### Nordsieck implementations

  - **`AN5`**: Fifth-order Adams method with fixed leading coefficient
  - **`JVODE`**: Variable order Adams/BDF method (experimental LSODE-style)
  - **`JVODE_Adams`**: JVODE configured for Adams methods
  - **`JVODE_BDF`**: JVODE configured for BDF methods

### Recommended alternatives

  - **For most applications**: Use `QNDF` or `FBDF` instead
  - **For stiff problems**: Prefer modern BDF implementations
  - **For research**: These methods are appropriate for experimental work

## Research and Development

These implementations serve as:

  - **Experimental testbed** for Nordsieck form algorithms
  - **Research platform** for discontinuity-aware methods
  - **Development basis** for future improved BDF implementations
  - **Educational examples** of alternative multistep representations

## Usage Guidelines

  - **Not recommended** for production applications
  - **Use FBDF or QNDF** for reliable multistep integration
  - **Consider these methods** only for research or experimental purposes
  - **Expect potentially lower performance** compared to modern alternatives

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqNordsieck", "AN5")
```

## Full list of solvers

```@docs
AN5
JVODE
JVODE_Adams
JVODE_BDF
```
