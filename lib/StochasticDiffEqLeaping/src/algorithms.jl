"""
    TauLeaping()

**TauLeaping: Basic Tau-Leaping Method (Jump-Diffusion)**

Basic tau-leaping method for approximating jump-diffusion processes by "leaping" over multiple potential jump events.

## Method Properties

  - **Problem type**: Jump-diffusion processes
  - **Approach**: Approximate multiple jumps per time step
  - **Time stepping**: Fixed tau approach
  - **Accuracy**: Depends on tau selection

## When to Use

  - Jump-diffusion systems with many small jumps
  - When exact jump simulation is computationally prohibitive
  - Chemical reaction networks with fast reactions
  - Population models with high birth-death rates
  - Initial exploration of jump-diffusion problems

## Algorithm Description

Approximates Poisson processes by assuming constant propensities over time interval tau, then sampling number of jumps from Poisson distribution.

## Tau Selection

Critical parameter: tau should be small enough that jump rates don't change significantly over [t, t+tau].

## References

  - Gillespie, D.T., "Approximate accelerated stochastic simulation of chemically reacting systems"
"""
@kwdef struct TauLeaping{QT} <: StochasticDiffEqJumpAdaptiveAlgorithm
    gamma::QT = 9 // 10
    qmax::QT = 10 // 1
end

"""
    CaoTauLeaping()

**CaoTauLeaping: Cao's Adaptive Tau-Leaping Method (Jump-Diffusion)**

Advanced tau-leaping method with adaptive tau selection and improved error control.

## Method Properties

  - **Problem type**: Jump-diffusion processes
  - **Approach**: Adaptive tau selection with error control
  - **Time stepping**: Adaptive tau based on error estimates
  - **Accuracy**: Superior to basic tau-leaping

## When to Use

  - Production jump-diffusion simulations requiring reliability
  - When adaptive tau selection is needed
  - Problems where basic TauLeaping gives poor accuracy
  - Chemical reaction networks requiring precise control

## Algorithm Features

  - Adaptive tau selection based on error estimates
  - Better stability and accuracy than basic tau-leaping
  - Automatic step size control
  - More sophisticated error estimation

## Tau Selection

Automatically adjusts tau based on:

  - Local error estimates
  - Jump rate variations
  - Solution stability requirements

## References

  - Cao, Y., Gillespie, D.T., Petzold, L.R., "Efficient step size selection for the tau-leaping method"
"""
@kwdef struct CaoTauLeaping{QT} <: StochasticDiffEqJumpAdaptiveAlgorithm
    gamma::QT = 9 // 10
end

"""
    ImplicitTauLeaping(; nlsolve=NLFunctional())

**ImplicitTauLeaping: First Order Implicit Tau-Leaping Method (Jump-Diffusion)**

An implicit (backward Euler) tau-leaping method for stiff chemical kinetic systems.
Uses backward Euler discretization to provide improved stability for systems with
fast reversible reactions or stiff rate constants.

## Method Properties

| Property              | Value       |
|:----------------------|:------------|
| Jacobian Required     | No          |
| Implicit              | Yes         |
| Adaptive              | No          |
| Stability             | A-stable    |
| Weak Order            | 1           |

## Mathematical Formulation

The method solves the implicit equation:

```math
X_{n+1} = X_n + ν ⋅ Poisson(dt ⋅ a(X_{n+1}))
```

which is approximated by:

```math
X_{n+1} = X_n + ν ⋅ k + dt ⋅ (drift(X_{n+1}) - drift(X_n))
```

where k ~ Poisson(dt * a(X_n)) and drift(u) = ν * a(u).

This corresponds to `ThetaTrapezoidalTauLeaping` with θ = 1 (fully implicit).

## Keyword Arguments

- `nlsolve`: Nonlinear solver algorithm (default: `NLFunctional()`).
  Options include `NLFunctional()`, `NLAnderson()`, and `NLNewton()`.

## Example

```julia
using StochasticDiffEq, JumpProcesses

# Define rate function and stoichiometry
rate(out, u, p, t) = (out[1] = 0.1*u[1]; out[2] = 0.05*u[2])
c(du, u, p, t, counts, mark) = (du[1] = -counts[1] + counts[2]; du[2] = counts[1] - counts[2])

rj = RegularJump(rate, c, 2)
prob = DiscreteProblem([100.0, 0.0], (0.0, 10.0))
jprob = JumpProblem(prob, Direct(), rj)

sol = solve(jprob, ImplicitTauLeaping(); dt=0.1)
```

## References

  - Rathinam, M., Petzold, L.R., Cao, Y., Gillespie, D.T., "Stiffness in stochastic
    chemically reacting systems: The implicit tau-leaping method", J. Chem. Phys.
    119, 12784 (2003)
"""
struct ImplicitTauLeaping{N} <: StochasticDiffEqJumpAdaptiveAlgorithm
    nlsolve::N
end
function ImplicitTauLeaping(; nlsolve = NLFunctional())
    return ImplicitTauLeaping(nlsolve)
end

"""
    ThetaTrapezoidalTauLeaping(; theta=0.5, max_iters=10, abstol=1e-8, reltol=1e-6)

**ThetaTrapezoidalTauLeaping: Implicit Weak Second Order Tau-Leaping Method (Jump-Diffusion)**

An implicit tau-leaping method achieving weak second order accuracy in the
large volume scaling. Uses fixed-point iteration to solve the implicit equation.
Based on the work of Hu, Li, and Min (2011) and Anderson and Mattingly (2011).

## Method Properties

  - **Problem type**: Jump-diffusion processes
  - **Order**: Weak order 2 (in the large volume scaling)
  - **Time stepping**: Fixed or adaptive tau
  - **Accuracy**: Superior to both Euler tau-leaping and midpoint tau-leaping
  - **Implicit treatment**: Uses nonlinear solver for implicit rate equations

## Parameters

  - `theta::Float64`: Implicitness parameter (default: 0.5)
    - Must be in range (0, 1)
    - theta = 0.5 gives trapezoidal method (recommended for balanced accuracy/stability)
    - theta = 1.0 gives backward Euler (maximum stability)
  - `max_iters::Int`: Maximum iterations for nonlinear solver (default: 10)
  - `abstol::Float64`: Absolute tolerance for convergence (default: 1e-8)
  - `reltol::Float64`: Relative tolerance for convergence (default: 1e-6)

## When to Use

  - When higher accuracy is needed compared to standard tau-leaping methods
  - Chemical reaction networks requiring weak second order accuracy
  - Systems where accurate mean and covariance estimates are important
  - When both Euler and midpoint tau-leaping provide insufficient accuracy
  - Stiff chemical systems where implicit treatment provides stability

## Algorithm Description

The method solves the implicit equation:

```math
X_{n+1} = X_n + ν⋅k + θ⋅dt⋅ν⋅(a(X_{n+1}) - a(X_n))
```

where `k ~ Poisson(dt⋅a(X_n))` are the jump counts.

This is solved using fixed-point iteration:
1. Generate Poisson jumps with rate `dt·a(X_n)`
2. Initialize `z = 0`
3. Iterate: `z_{new} = θ·dt·(drift(X_n + ν·k + z) - drift(X_n))`
4. Final state: `X_{n+1} = X_n + ν·k + z`

## Convergence Properties

The local truncation error for covariance is O(τ³V⁻¹) when τ = V^(-β) for 0 < β < 1
and system size V → ∞, which is higher order than both Euler and midpoint methods.

## References

  - Hu, Y., Li, T., Min, B., "A weak second order tau-leaping method for chemical
    kinetic systems", J. Chem. Phys. 135, 024113 (2011)
  - Anderson, D.F., Mattingly, J.C., "A weak trapezoidal method for a class of
    stochastic differential equations", Comm. Math. Sci. 9, 301 (2011)
"""
struct ThetaTrapezoidalTauLeaping{T, N} <: StochasticDiffEqJumpAdaptiveAlgorithm
    theta::T
    nlsolve::N
end
function ThetaTrapezoidalTauLeaping(; theta = 0.5, nlsolve = NLFunctional())
    return ThetaTrapezoidalTauLeaping(theta, nlsolve)
end
