```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqRKN

Runge-Kutta-Nyström (RKN) methods for solving second-order differential equations of the form `d²u/dt² = f(u, du/dt, t)`. These methods are specifically designed for second-order ODEs and can be more efficient than converting to first-order systems when the problem has this natural structure.

## Key Properties

RKN methods provide:

  - **Direct integration** of second-order ODEs without conversion to first-order systems
  - **Efficient handling** of problems where velocity dependence is minimal or absent
  - **Specialized variants** for different types of second-order problems
  - **High-order accuracy** with methods up to 12th order available
  - **Trigonometrically-fitted variants** for oscillatory problems
  - **Improved efficiency** for second-order problems compared to first-order conversions

## When to Use RKN Methods

These methods are recommended for:

  - **Second-order differential equations** with natural `d²u/dt² = f(u, du/dt, t)` structure
  - **Classical mechanics problems** (Newton's equations of motion)
  - **Oscillatory second-order systems** (springs, pendulums, wave equations)
  - **Problems where velocity dependence is weak** or absent
  - **Celestial mechanics** and orbital dynamics
  - **Vibration analysis** and structural dynamics
  - **Wave propagation** problems in their natural second-order form

## Problem Types

### Velocity-independent problems: `d²u/dt² = f(u, t)`

When the acceleration depends only on position (and possibly time):

  - **More efficient specialized methods** available
  - **Classical examples**: gravitational problems, springs with `F = -kx`

### Velocity-dependent problems: `d²u/dt² = f(u, du/dt, t)`

When acceleration depends on both position and velocity:

  - **General RKN methods** handle the full dependence
  - **Examples**: damped oscillators, air resistance problems

## Solver Selection Guide

### General-purpose RKN methods

  - **`Nystrom4`**: Fourth-order method for general second-order ODEs
  - **`IRKN3`**, **`IRKN4`**: Improved Runge-Kutta-Nyström methods

### Velocity-independent specialized methods

  - **`Nystrom4VelocityIndependent`**: Fourth-order for `d²u/dt² = f(u, t)`
  - **`Nystrom5VelocityIndependent`**: Fifth-order for `d²u/dt² = f(u, t)`

### Velocity-dependent methods

  - **`FineRKN4`**: Fourth-order allowing velocity dependence
  - **`FineRKN5`**: Fifth-order allowing velocity dependence

### High-order Dormand-Prince variants

  - **`DPRKN4`**: Fourth-order Dormand-Prince RKN
  - **`DPRKN5`**: Fifth-order Dormand-Prince RKN
  - **`DPRKN6`**: Sixth-order with free interpolant
  - **`DPRKN6FM`**: Sixth-order with optimized error coefficients
  - **`DPRKN8`**: Eighth-order for high accuracy
  - **`DPRKN12`**: Twelfth-order for extreme precision

### Trigonometrically-fitted methods for oscillatory problems

  - **`ERKN4`**: Embedded 4(3) pair for periodic problems
  - **`ERKN5`**: Embedded 5(4) pair for periodic problems
  - **`ERKN7`**: Higher-order embedded pair for oscillatory systems

### Specialized methods

  - **`RKN4`**: Fourth-order for linear inhomogeneous second-order IVPs

## Advantages Over First-Order Conversion

  - **Natural problem structure**: Preserves the physical meaning of the equations
  - **Specialized optimizations**: Methods can exploit second-order structure

To be able to access the solvers in `OrdinaryDiffEqRKN`, you must first install them use the Julia package manager:

```julia
using Pkg
Pkg.add("OrdinaryDiffEqRKN")
```

This will only install the solvers listed at the bottom of this page.
If you want to explore other solvers for your problem,
you will need to install some of the other libraries listed in the navigation bar on the left.

## Example usage

```julia
using OrdinaryDiffEqRKN
function HH_acceleration!(dv, v, u, p, t)
    x, y = u
    dx, dy = dv
    dv[1] = -x - 2 * x * y
    dv[2] = y^2 - y - x^2
end
initial_positions = [0.0, 0.1]
initial_velocities = [0.5, 0.0]
tspan = (0.0, 1.0)
prob = SecondOrderODEProblem(HH_acceleration!, initial_velocities, initial_positions, tspan)
sol = solve(prob, Nystrom4(), dt = 1 / 10)
```

## Full list of solvers

```@docs
IRKN3
IRKN4
Nystrom4
Nystrom4VelocityIndependent
Nystrom5VelocityIndependent
FineRKN4
FineRKN5
DPRKN4
DPRKN5
DPRKN6
DPRKN6FM
DPRKN8
DPRKN12
ERKN4
ERKN5
ERKN7
RKN4
```
