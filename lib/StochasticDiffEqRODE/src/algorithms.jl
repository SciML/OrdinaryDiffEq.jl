"""
    RandomEM()

**RandomEM: Random Euler Method (RODE)**

Euler method for Random Ordinary Differential Equations (RODEs) with random parameters.

## Method Properties

  - **Problem type**: Random ODEs (RODEs)
  - **Strong Order**: 1.0 (for deterministic part)
  - **Randomness**: Handles random parameters, not Brownian motion
  - **Time stepping**: Fixed step size

## When to Use

  - Random ODEs with random parameters but no Brownian motion
  - Uncertainty quantification with parameter randomness
  - Problems with random coefficients or initial conditions
  - Monte Carlo simulation of deterministic systems with random inputs

## RODE vs SDE

  - **RODE**: Random parameters, deterministic evolution
  - **SDE**: Fixed parameters, stochastic (Brownian) evolution

## References

  - Random ordinary differential equation methods
"""
struct RandomEM <: StochasticDiffEqRODEAlgorithm end

"""
    RandomHeun()

**RandomHeun: Random Heun Method (RODE)**

Heun method for Random Ordinary Differential Equations with improved accuracy.

## Method Properties

  - **Problem type**: Random ODEs (RODEs)
  - **Strong Order**: 2.0 (for deterministic part)
  - **Randomness**: Handles random parameters
  - **Time stepping**: Fixed step size

## When to Use

  - RODEs requiring higher accuracy than RandomEM
  - When computational cost per step is acceptable
  - Random parameter problems needing second-order accuracy

## References

  - Higher-order methods for random ODEs
"""
struct RandomHeun <: StochasticDiffEqRODEAlgorithm end

"""
    RandomTamedEM()

**RandomTamedEM: Tamed Random Euler Method (RODE)**

Tamed Euler method for RODEs with potentially explosive behavior.

## Method Properties

  - **Problem type**: Random ODEs with potential blow-up
  - **Approach**: Taming to prevent numerical explosion
  - **Stability**: Enhanced stability for unstable random systems
  - **Time stepping**: Fixed step size with taming

## When to Use

  - RODEs that may exhibit explosive growth
  - When RandomEM gives unstable or explosive solutions
  - Random systems with strong nonlinearities
  - Problems requiring enhanced numerical stability

## Taming Mechanism

Applies taming technique to prevent numerical blow-up while maintaining accuracy for well-behaved solutions.

## References

  - Tamed methods for random differential equations
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
