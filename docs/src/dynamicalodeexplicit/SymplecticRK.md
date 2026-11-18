```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqSymplecticRK

Symplectic integrators are specialized methods for solving Hamiltonian systems and second-order differential equations that preserve important geometric properties of the phase space. These methods are essential for long-time integration of conservative mechanical systems.

## Key Properties

Symplectic integrators provide:

  - **Exact conservation of symplectic structure** in phase space
  - **Bounded energy error** over long time periods
  - **Excellent long-time stability** without secular drift
  - **Preservation of periodic orbits** and other geometric structures
  - **Linear energy drift** instead of quadratic (much better than standard methods)

## When to Use Symplectic Methods

Symplectic integrators are essential for:

  - **Hamiltonian systems** and conservative mechanical problems
  - **Molecular dynamics** and N-body simulations
  - **Celestial mechanics** and orbital computations
  - **Plasma physics** and charged particle dynamics
  - **Long-time integration** where energy conservation is critical
  - **Oscillatory problems** requiring preservation of periodic structure
  - **Classical mechanics problems** with known analytical properties

## Mathematical Background

For a Hamiltonian system with energy `H(p,q)`, symplectic integrators preserve the symplectic structure `dp ∧ dq`. While standard integrators have energy error growing quadratically over time, symplectic methods maintain bounded energy with only linear drift, making them superior for long-time integration.

## Solver Selection Guide

### First-order methods

  - **`SymplecticEuler`**: First-order, simplest symplectic method. Only recommended when the dynamics function `f` is not differentiable.

### Second-order methods

  - **`McAte2`**: Optimized second-order McLachlan-Atela method, **recommended for most applications**
  - **`VelocityVerlet`**: Second-order, common choice for molecular dynamics but less efficient in terms of accuracy than McAte2
  - **`VerletLeapfrog`**: Second-order, kick-drift-kick formulation
  - **`LeapfrogDriftKickDrift`**: Alternative second-order leapfrog
  - **`PseudoVerletLeapfrog`**: Modified Verlet scheme

### Third-order methods

  - **`Ruth3`**: Third-order method
  - **`McAte3`**: Optimized third-order McLachlan-Atela method

### Fourth-order methods

  - **`CandyRoz4`**: Fourth-order method
  - **`McAte4`**: Fourth-order McLachlan-Atela (requires quadratic kinetic energy)
  - **`CalvoSanz4`**: Optimized fourth-order method (requires quadratic kinetic energy)
  - **`McAte42`**: Alternative fourth-order method (requires quadratic kinetic energy)

### Higher-order methods

  - **`McAte5`**: Fifth-order McLachlan-Atela method
  - **`Yoshida6`**: Sixth-order method
  - **`KahanLi6`**: Optimized sixth-order method
  - **`McAte8`**: Eighth-order McLachlan-Atela method
  - **`KahanLi8`**: Optimized eighth-order method
  - **`SofSpa10`**: Tenth-order method for highest precision

## Method Selection Guidelines

  - **For most applications**: `McAte2` (second-order, optimal efficiency)
  - **For molecular dynamics (common choice)**: `VelocityVerlet` (less efficient than McAte2 but widely used)
  - **For non-differentiable dynamics**: `SymplecticEuler` (first-order, only when necessary)
  - **For computational efficiency**: `McAte2` or `McAte3`

### Stage times

The tableau methods (`PseudoVerletLeapfrog`, `McAte2`, `Ruth3`, `McAte3`, `CandyRoz4`,
`McAte4`, `CalvoSanz4`, `McAte42`, `McAte5`, `Yoshida6`, `KahanLi6`, `McAte8`, `KahanLi8`,
`SofSpa10`) evaluate kick `i` of a step from `t` to `t + dt` at the time the position has
reached, `t + (b[1] + ... + b[i]) * dt`, and drift `i` at the time the velocity has reached,
`t + (a[1] + ... + a[i-1]) * dt`; the last evaluations of a step are at `t + dt`. This keeps
the order of the method for time-dependent forces. A time-dependent drift (`f2` depending
on `t`) is outside the quadratic-kinetic-energy assumption of `McAte4`, `CalvoSanz4` and
`McAte5`, which then drop to third order; the package tests check the stated order with a
time-dependent drift for `PseudoVerletLeapfrog`, `McAte2`, `Ruth3`, `McAte3`, `CandyRoz4`
and `McAte42`.

Methods with negative coefficients evaluate some stages outside `[t, t + dt]`, so on the
first or last step outside the integration interval:

  - the force `f1` (kick) for `Ruth3` (up to `t + 25/24 dt`) and `CandyRoz4` (from
    about `t - 0.351207 dt` to `t + 1.351207 dt`);
  - the drift `f2` for `McAte4`, `CalvoSanz4`, `McAte42`, `McAte5`, `Yoshida6` and
    `McAte8`, which matters only when `f2` depends on `t`.

If the right-hand side is defined only on the integration interval, extend it slightly
beyond the interval or use a method whose stage times all stay inside `[t, t + dt]`:
`PseudoVerletLeapfrog`, `McAte2`, `McAte3`, `KahanLi6`, `KahanLi8` or `SofSpa10`.

### Important Note on Chaotic Systems

Most N-body problems (molecular dynamics, astrophysics) are chaotic systems where solutions diverge onto shadow trajectories. In such cases, **higher-order methods provide no practical advantage** because the true error remains O(1) for sufficiently long integrations - exactly the scenarios where symplectic methods are most needed. The geometric properties preserved by symplectic integrators are more important than high-order accuracy for chaotic systems.

For more information on chaos and accuracy in numerical integration, see: [How Chaotic is Chaos? How Some AI for Science (SciML) Papers are Overstating Accuracy Claims](https://www.stochasticlifestyle.com/how-chaotic-is-chaos-how-some-ai-for-science-sciml-papers-are-overstating-accuracy-claims/)

## Installation

To be able to access the solvers in `OrdinaryDiffEqSymplecticRK`, you must first install them use the Julia package manager:

```julia
using Pkg
Pkg.add("OrdinaryDiffEqSymplecticRK")
```

This will only install the solvers listed at the bottom of this page.
If you want to explore other solvers for your problem,
you will need to install some of the other libraries listed in the navigation bar on the left.

## Example usage

```julia
using OrdinaryDiffEqSymplecticRK
function HH_acceleration!(dv, v, u, p, t)
    x, y = u
    dv[1] = dx = -x - 2 * x * y
    dv[2] = dy = y^2 - y - x^2
    return
end
initial_positions = [0.0, 0.1]
initial_velocities = [0.5, 0.0]
tspan = (0.0, 1.0)
prob = SecondOrderODEProblem(HH_acceleration!, initial_velocities, initial_positions, tspan)
sol = solve(prob, KahanLi8(), dt = 1 / 10)
```

## Full list of solvers

```@docs
SymplecticEuler
VelocityVerlet
VerletLeapfrog
LeapfrogDriftKickDrift
PseudoVerletLeapfrog
McAte2
Ruth3
McAte3
CandyRoz4
McAte4
CalvoSanz4
McAte42
McAte5
Yoshida6
KahanLi6
McAte8
KahanLi8
SofSpa10
```
