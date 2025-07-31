```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqLowOrderRK

Low-order explicit Runge-Kutta methods for non-stiff differential equations. **Most of the time, you should use [`Tsit5`](@ref OrdinaryDiffEqTsit5)**, which is the most common and efficient low-order RK method. The alternative methods provided here are for special circumstances where Tsit5 is not suitable.

## Key Properties

Low-order explicit RK methods offer:

  - **Computational efficiency** at higher tolerances (>1e-6)
  - **Robust error control** for difficult non-stiff problems
  - **Specialized interpolation properties** for applications requiring dense output
  - **Lower-order derivatives** requirements for non-smooth functions
  - **Good performance** for specific problem types

## When to Use Alternative Low-Order RK Methods

Choose these methods instead of Tsit5 when:

  - **ODE function `f` is not differentiable to 5th order** - use lower-order methods (the more discontinuous, the lower the order needed)
  - **Heavy use of interpolations** - OwrenZen methods have superior interpolation convergence
  - **Delay differential equations** - OwrenZen methods are most efficient (see SciMLBenchmarks)
  - **Very high tolerances** (>1e-3) - BS3 is more efficient than Tsit5
  - **Quadratic polynomial ODEs** - SIR54 is optimized for these systems
  - **Educational purposes** - simpler methods for understanding algorithms

## Solver Selection Guide

### Primary recommendation

**For most problems, use [`Tsit5`](@ref OrdinaryDiffEqTsit5) instead of these methods.**

### High tolerances (>1e-3)

  - **`BS3`**: Third-order Bogacki-Shampine method, most efficient for very high tolerances

### Superior interpolation needs

  - **`OwrenZen3`**: Third-order with excellent interpolation convergence
  - **`OwrenZen5`**: Fifth-order with excellent interpolation, **optimal for DDEs**
  - **`OwrenZen4`**: Fourth-order interpolation-optimized method

### Non-smooth or discontinuous ODEs

  - **`BS3`**: Third-order for mildly non-smooth functions
  - **`Heun`**: Second-order for more discontinuous functions (not generally recommended)
  - **`Euler`**: First-order for highly discontinuous problems

### Robust error control alternatives

  - **`BS5`**: Fifth-order with very robust error estimation
  - **`DP5`**: Fifth-order Dormand-Prince method, classical alternative to Tsit5

### Specialized applications

  - **`RK4`**: Fourth-order with special residual error control, good for DDEs. **Note**: Uses adaptive timestepping by default - set `adaptive=false` in `solve()` for traditional fixed-step RK4
  - **`SIR54`**: Fifth-order optimized for ODEs defined by quadratic polynomials (e.g., SIR-type epidemiological models)
  - **`Stepanov5`**: Fifth-order method with enhanced stability properties and optimized error constants
  - **`Ralston`**: Second-order with optimized error constants

### Periodic and oscillatory problems

  - **`Anas5`**: Fifth-order optimized for periodic problems with minimal phase error
  - **`FRK65`**: Sixth-order zero dissipation method for oscillatory problems

### Advanced specialized methods

  - **`RKO65`**: Sixth-order optimized method
  - **`MSRK5`**, **`MSRK6`**: Multi-stage methods for specific applications
  - **`PSRK4p7q6`**, **`PSRK3p5q4`**, **`PSRK3p6q5`**: Pseudo-symplectic methods
  - **`Alshina2`**, **`Alshina3`**, **`Alshina6`**: Methods with optimized parameters

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqLowOrderRK", "BS3")
```

## Full list of solvers

```@docs
Euler
Heun
Ralston
Midpoint
RK4
BS3
OwrenZen3
OwrenZen4
OwrenZen5
BS5
DP5
Anas5
RKO65
FRK65
RKM
MSRK5
MSRK6
PSRK4p7q6
PSRK3p5q4
PSRK3p6q5
Stepanov5
SIR54
Alshina2
Alshina3
Alshina6
```
