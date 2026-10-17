"""
    RandomEM()

**RandomEM: Random Euler Method (RODE)**

Euler method for Random Ordinary Differential Equations. Each step advances `u` by
`dt * f(u, p, t, W)` with the driving noise `W` taken at the start of the step. Any
noise process can drive the problem; a `RODEProblem` given no noise is driven by a
Wiener process.

## Method Properties

  - **Problem type**: Random ODEs (RODEs)
  - **Strong order**: 1 on Wiener and other semimartingale noise, under the conditions
    of Kloeden and Rosa
  - **Time stepping**: Fixed step size
  - **Right-hand side evaluations**: 1 per step

## When to Use

The default RODE method. Higher-order classical tableaus do not raise the rate on a
Wiener-driven RODE, so this is the reference for cost; `RandomTaylor15` reaches past
order 1 when the path is stored on a grid finer than the solver steps.

## References

  - Kloeden and Rosa, Strong order-one convergence of the Euler method for random
    ordinary differential equations driven by semi-martingale noises, ESAIM: M2AN 59
    (2025), arXiv:2306.15418.
"""
struct RandomEM <: StochasticDiffEqRODEAlgorithm end

"""
    RandomHeun()

**RandomHeun: Random Heun Method (RODE)**

Two-stage Heun method for Random Ordinary Differential Equations. The predictor uses
`W` at the start of the step and the corrector uses `W` at the end. When `f` does not
depend on `W` this is the classical second-order Heun method.

## Method Properties

  - **Problem type**: Random ODEs (RODEs)
  - **Strong order**: 1 measured on Wiener noise; 2 when `f` does not depend on the
    noise
  - **Time stepping**: Fixed step size
  - **Right-hand side evaluations**: 2 per step

## When to Use

  - When `f` depends weakly on the noise, so the second stage lowers the error constant
  - On a Wiener-driven RODE it measures the same order 1 as `RandomEM` at twice the cost
    per step, so it does not raise the rate there
"""
struct RandomHeun <: StochasticDiffEqRODEAlgorithm end

"""
    RandomTamedEM()

**RandomTamedEM: Tamed Random Euler Method (RODE)**

Tamed Euler method for Random Ordinary Differential Equations. Each step advances `u`
by `dt * k / (1 + dt * norm(k))` with `k = f(u, p, t, W)` and `W` taken at the start of
the step, so no step moves `u` by more than about 1 in the Euclidean norm however large
`f` is. The bound is in the units of `u`, so rescaling `u` changes the numerical
solution.

## Method Properties

  - **Problem type**: Random ODEs (RODEs)
  - **Strong order**: 1 measured on Wiener noise
  - **Time stepping**: Fixed step size
  - **Right-hand side evaluations**: 1 per step

## When to Use

  - RODEs whose `f` grows superlinearly, where `RandomEM` can blow up

## References

  - Hutzenthaler, Jentzen and Kloeden, Strong convergence of an explicit numerical
    method for SDEs with non-globally Lipschitz continuous coefficients, Ann. Appl.
    Probab. 22 (2012), DOI 10.1214/11-AAP803. The taming factor comes from there; its
    convergence result is for SDEs, not RODEs.
"""
struct RandomTamedEM <: StochasticDiffEqRODEAlgorithm end

"""
    RandomTaylor15()

**RandomTaylor15: Derivative-free order 1.5 Taylor method (RODE)**

Order 1.5 scheme for Random Ordinary Differential Equations driven by a Wiener process
supplied as a stored path. The derivatives of the Taylor scheme are replaced by finite
differences of the right-hand side, so only evaluations of `f` are needed. The drift is
advanced by a Heun step with the noise held at its value at `t`, and the noise enters
through integrals of the path over the step, taken from the supplied path rather than
from their Brownian expectations. The step therefore reduces to Heun's method when `f`
does not depend on `W`, and its order does not depend on the amplitude of the path.

## Method Properties

  - **Problem type**: RODEs driven by a `NoiseGrid` with scalar values
  - **Pathwise order**: 1.5
  - **Time stepping**: Fixed step size
  - **Right-hand side evaluations**: 4 per step
  - **Noise usage**: reads the driving path between the step endpoints

## When to Use

The step uses integrals of the driving path over `[t, t+dt]`, so the path must be
resolved more finely than the solver steps. That happens when the noise is measured
data or is generated on a fine grid and the solver is stepped coarsely. With a path
that is only known at the solver's own steps there is no sub-step information to use
and `RandomEM` is the appropriate method; the integrals then collapse to the endpoint
rule and the order drops to 1, which is warned about when the cache is built.

The finite differences perturb the noise argument by `sqrt(dt)`, so the error constant
carries the third derivative of `f` in `W`. For a right-hand side oscillating in `W` at
frequency `a` the asymptotic rate is reached once `a^2 * dt` is below about 1, and the
measured rate is lower on coarser steps.

## References

  - Asai, Numerical Methods for Random Ordinary Differential Equations and their
    Applications in Biology and Medicine, PhD thesis, Goethe University Frankfurt, 2016,
    equation (3.24), with the step integrals taken from the path and the drift advanced
    by a Heun step rather than through the second difference.
"""
struct RandomTaylor15 <: StochasticDiffEqRODEAlgorithm end

"""
    BAOAB(; gamma = 1.0, scale_noise = true)

**BAOAB: Langevin Dynamics Integrator (Specialized)**

Specialized integrator for Langevin dynamics in molecular dynamics simulations, particularly effective for configurational sampling.

## Method Properties

  - **Problem type**: Langevin dynamics (second-order SDEs)
  - **Structure**: Position-velocity formulation
  - **Sampling**: Designed for equilibrium sampling
  - **Time stepping**: Fixed step size
  - **Conservation**: Preserves equilibrium distributions

## Parameters

  - `gamma::Real = 1.0`: Friction coefficient
  - `scale_noise::Bool = true`: Whether to scale noise appropriately

## System Structure

Designed for Langevin systems:

```math
\\begin{align*}
du &= v \\, dt \\\\
dv &= f(v,u) \\, dt - γv \\, dt + g(u) \\sqrt{2γ} \\, dW
\\end{align*}
```

where:

  - ``u``: position coordinates
  - ``v``: velocity coordinates
  - ``γ``: friction coefficient
  - ``f(v,u)``: force function
  - ``g(u)``: noise scaling function

## When to Use

  - Molecular dynamics simulations with Langevin thermostat
  - Configurational sampling of molecular systems
  - Equilibrium sampling from canonical ensemble
  - Second-order SDEs with damping and noise

## Algorithm Features

  - BAOAB splitting: B(kick) - A(drift) - O(Ornstein-Uhlenbeck) - A(drift) - B(kick)
  - Preserves correct equilibrium distribution
  - Robust and efficient for molecular sampling
  - Well-suited for long-time integration

## References

  - Leimkuhler B., Matthews C., "Robust and efficient configurational molecular sampling via Langevin dynamics", J. Chem. Phys. 138, 174102 (2013)
"""
struct BAOAB{T} <: StochasticDiffEqAlgorithm
    gamma::T
    scale_noise::Bool
end
BAOAB(; gamma = 1.0, scale_noise = true) = BAOAB(gamma, scale_noise)
