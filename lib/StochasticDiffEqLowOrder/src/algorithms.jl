"""
    EM(; split = true)

Select the fixed-step Euler-Maruyama method for an Itô SDE.

# Keywords

- `split = true`: Whether to use split-step handling.

# Examples

```julia
sol = solve(prob, EM(); dt = 0.01)
```

Use `split = false` when the problem should be advanced without split-step
handling.
"""
struct EM{split} <: StochasticDiffEqAlgorithm end
EM(split = true) = EM{split}()

"""
    SplitEM()

Select the fixed-step split Euler-Maruyama method for an Itô SDE with a
splittable drift.

# Examples

```julia
sol = solve(prob, SplitEM(); dt = 0.01)
```
"""
struct SplitEM <: StochasticDiffEqAlgorithm end

"""
    EulerHeun()

**EulerHeun: Stochastic Euler-Heun Method**

A two-stage predictor-corrector (Heun-type) method for Stratonovich SDEs, the Stratonovich
analogue of Euler-Maruyama. An Euler step forms the predictor `ũ`, then the drift and
diffusion are re-evaluated there and trapezoidally averaged to form the update:

  - predictor: `ũ = uₙ + f(uₙ)·Δt + g(uₙ)·ΔW`
  - corrector: `uₙ₊₁ = uₙ + ½(f(uₙ) + f(ũ))·Δt + ½(g(uₙ) + g(ũ))·ΔW`

This corresponds to the improved-Euler scheme of Roberts (2012). Note this is a genuine
two-stage scheme (two drift and two diffusion evaluations per step), not a single-stage
predictor-corrector.

## Method Properties

  - **Strong Order**: 1.0 for Stratonovich SDEs with commutative noise (e.g. scalar,
    diagonal, or additive); 1/2 for general non-commutative noise, which is the value
    reported by `alg_order` (the same convention as `EM`)
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed step size
  - **Noise types**: General (scalar, diagonal, non-diagonal)
  - **SDE interpretation**: Stratonovich

## When to Use

  - For Stratonovich SDEs where a simple, low-cost method is sufficient
  - As the Stratonovich counterpart to `EM` for Itô problems
  - When adaptive time stepping is not required; see `LambaEulerHeun` for an adaptive variant

## References

  - Roberts, A.J., "Modify the improved Euler scheme to integrate stochastic differential
    equations", arXiv:1210.0933 (2012)
  - Kloeden, P.E., Platen, E., Numerical Solution of Stochastic Differential Equations,
    Springer, Berlin Heidelberg, p. 373 (1992)
"""
struct EulerHeun <: StochasticDiffEqAlgorithm end

"""
    LambaEM(; split = true)

Select the adaptive Lamba Euler-Maruyama method for an Itô SDE.

# Keywords

- `split = true`: Whether to use split-step handling.

# Examples

```julia
sol = solve(prob, LambaEM(); reltol = 1.0e-3, abstol = 1.0e-3)
```
"""
struct LambaEM{split} <: StochasticDiffEqAdaptiveAlgorithm end
LambaEM(split = true) = LambaEM{split}()

"""
    LambaEulerHeun()

Select the adaptive Lamba Euler-Heun method for a Stratonovich SDE.

# Examples

```julia
sol = solve(prob, LambaEulerHeun(); reltol = 1.0e-3, abstol = 1.0e-3)
```
"""
struct LambaEulerHeun <: StochasticDiffEqAdaptiveAlgorithm end

"""
    SimplifiedEM()

Select the simplified fixed-step Euler-Maruyama method for weak SDE
approximations.

# Examples

```julia
sol = solve(prob, SimplifiedEM(); dt = 0.01)
```
"""
struct SimplifiedEM <: StochasticDiffEqAlgorithm end

"""
    RKMil(; interpretation = SciMLBase.AlgorithmInterpretation.Ito)

Select the adaptive Runge-Kutta Milstein method for scalar or diagonal-noise
SDEs.

# Keywords

- `interpretation = SciMLBase.AlgorithmInterpretation.Ito`: Stochastic
  interpretation used by the method.

# Examples

```julia
sol = solve(prob, RKMil(); reltol = 1.0e-3, abstol = 1.0e-3)
```
"""
struct RKMil{interpretation} <: StochasticDiffEqAdaptiveAlgorithm end
RKMil(; interpretation = SciMLBase.AlgorithmInterpretation.Ito) = RKMil{interpretation}()

"""
    RKMilCommute(;
        interpretation = SciMLBase.AlgorithmInterpretation.Ito,
        ii_approx = IICommutative()
    )

Select the adaptive Runge-Kutta Milstein method for commutative-noise SDEs.

# Fields

- `interpretation`: Stochastic interpretation used by the method.
- `ii_approx`: Iterated-integral approximation for commutative noise.

# Keywords

- `interpretation = SciMLBase.AlgorithmInterpretation.Ito`: Stochastic
  interpretation used by the method.
- `ii_approx = IICommutative()`: Iterated-integral approximation for
  commutative noise.

# Examples

```julia
sol = solve(prob, RKMilCommute(); reltol = 1.0e-3, abstol = 1.0e-3)
```
"""
struct RKMilCommute{T} <: StochasticDiffEqAdaptiveAlgorithm
    interpretation::SciMLBase.AlgorithmInterpretation.T
    ii_approx::T
end
function RKMilCommute(; interpretation = SciMLBase.AlgorithmInterpretation.Ito, ii_approx = IICommutative())
    return RKMilCommute(interpretation, ii_approx)
end

"""
    PCEuler(ggprime; theta = 1 / 2, eta = 1 / 2)

Select a fixed-step predictor-corrector Euler method for an Itô SDE, supplying
the derivative of the diffusion coefficient as `ggprime`.

# Arguments

- `ggprime`: Derivative of the diffusion coefficient.

# Keywords

- `theta = 1 / 2`: Drift evaluation weight.
- `eta = 1 / 2`: Diffusion evaluation weight.

# Examples

```julia
sol = solve(prob, PCEuler(ggprime); dt = 0.01)
```
"""
struct PCEuler{T <: Real, F} <: StochasticDiffEqAlgorithm
    theta::T
    eta::T
    ggprime::F
end
PCEuler(ggprime; theta = 1 / 2, eta = 1 / 2) = PCEuler(theta, eta, ggprime)
