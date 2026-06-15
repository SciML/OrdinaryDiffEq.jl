struct EM{split} <: StochasticDiffEqAlgorithm end
EM(split = true) = EM{split}()

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

struct LambaEM{split} <: StochasticDiffEqAdaptiveAlgorithm end
LambaEM(split = true) = LambaEM{split}()

struct LambaEulerHeun <: StochasticDiffEqAdaptiveAlgorithm end

struct SimplifiedEM <: StochasticDiffEqAlgorithm end

struct RKMil{interpretation} <: StochasticDiffEqAdaptiveAlgorithm end
RKMil(; interpretation = SciMLBase.AlgorithmInterpretation.Ito) = RKMil{interpretation}()

struct RKMilCommute{T} <: StochasticDiffEqAdaptiveAlgorithm
    interpretation::SciMLBase.AlgorithmInterpretation.T
    ii_approx::T
end
function RKMilCommute(; interpretation = SciMLBase.AlgorithmInterpretation.Ito, ii_approx = IICommutative())
    return RKMilCommute(interpretation, ii_approx)
end

struct PCEuler{T <: Real, F} <: StochasticDiffEqAlgorithm
    theta::T
    eta::T
    ggprime::F
end
PCEuler(ggprime; theta = 1 / 2, eta = 1 / 2) = PCEuler(theta, eta, ggprime)
