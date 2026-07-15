```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqNewmark

Newmark-β and generalized-α methods for second-order ODEs, typically in mass-matrix form. These are second order time integrators, advancing the displacement and velocity with Newmark updates and solve a nonlinear residual for the acceleration directly. They arose from time stepping in structural and computational mechanics and are designed for such systems.

## Key Properties

These methods provide:

  - **Direct integration** of second-order ODEs, including mass-matrix form
  - **Unconditionally stable** parameter choices for structural dynamics
  - **Controllable high-frequency damping** via generalized-α's (`ρ∞`)
  - **Newmark-β and HHT-α** as special cases of generalized-α
  - **Adaptive time stepping** with the Zienkiewicz–Xie local truncation error estimate (Newmark-β)

## When to Use Newmark Methods

These methods are recommended for:

  - **Structural dynamics** and finite-element/Galerkin vibration problems
  - **Second-order ODEs** of the form `M a = f(u, v, t)`
  - **Problems needing algorithmic damping** of spurious high-frequency modes
  - **Mass-matrix second-order systems** arising from spatial discretizations
  - **Cases where RKN/symplectic methods are not the right fit** (implicit structural integrators with dissipation control)

## Mathematical Background

The generalized-α method evaluates the equations of motion at interpolated states:

`M * aₙ₊αₘ = f(uₙ₊αf, vₙ₊αf, tₙ₊αf)`

with

  - `aₙ₊αₘ = (1 - αₘ) * aₙ₊₁ + αₘ * aₙ`
  - `uₙ₊αf = (1 - αf) * uₙ₊₁ + αf * uₙ`
  - `vₙ₊αf = (1 - αf) * vₙ₊₁ + αf * vₙ`

and the standard Newmark updates for `uₙ₊₁` and `vₙ₊₁`. Setting `αₘ = αf = 0` recovers Newmark-β; setting `αₘ = 0` recovers HHT-α.

## Solver Selection Guide

  - **`NewmarkBeta`**: Classical Newmark-β. Default `β = 1/4`, `γ = 1/2` (average acceleration, second-order when `γ = 1/2`).
  - **`GeneralizedAlpha`**: Preferred when controllable high-frequency damping is needed.
      - **`GeneralizedAlpha(; rho_inf)`**: Recommended parameterization. `ρ∞ = 1` gives no algorithmic damping (equivalent to undamped Newmark); `ρ∞ = 0` gives maximum damping. Always sets unconditionally stable in [0, 1].
      - **`GeneralizedAlpha(; alpha_hht)`**: HHT-α convenience (`α ∈ [-1/3, 0]`). Also unconditionally stable in that range.
      - **`GeneralizedAlpha(αm, αf, β, γ)`**: Explicit four-parameter construction. Can break unconditionally stablity if parameters are chosen poorly.

### Unconditional stability

**Newmark-β** is unconditionally stable when `γ ≥ 1/2` and `β ≥ (1/4) * (γ + 1/2)^2`.

The default `β = 1/4`, `γ = 1/2` (average acceleration) sits on this bound and is second-order. `γ = 1/2` with `β < 1/4` (e.g. central difference `β = 0`) is only conditionally stable. `γ > 1/2` adds numerical damping but drops the method to first order.

**Generalized-α** is unconditionally stable when `αₘ ≤ αf ≤ 1/2` and `β ≥ (1/4) * (1/2 + αf - αₘ)^2`.

Second-order accuracy further requires `γ = 1/2 - αₘ + αf`. The `rho_inf` and `alpha_hht` constructors enforce these choices; the four-parameter form asserts the stability inequalities above at construction, so values that violate them will error rather than silently run unstably.

## Installation

To be able to access the solvers in `OrdinaryDiffEqNewmark`, you must first install them use the Julia package manager:

```julia
using Pkg
Pkg.add("OrdinaryDiffEqNewmark")
```

This will only install the solvers listed at the bottom of this page.
If you want to explore other solvers for your problem,
you will need to install some of the other libraries listed in the navigation bar on the left.

## Example usage

```julia
using OrdinaryDiffEqNewmark
function f1!(dv, v, u, p, t)
    dv .= -u
end
function f2!(du, v, u, p, t)
    du .= v
end
v0 = ones(2)
u0 = zeros(2)
prob = DynamicalODEProblem(f1!, f2!, v0, u0, (0.0, 5.0))
sol = solve(prob, NewmarkBeta(), dt = 0.1)
sol_ga = solve(prob, GeneralizedAlpha(; rho_inf = 0.8), dt = 0.1)
```

## Full list of solvers

```@docs
NewmarkBeta
GeneralizedAlpha
```
