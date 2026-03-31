"""
    RKMilGeneral(; interpretation=Ito, ii_approx=IILevyArea(), c=1, p=nothing, dt=nothing)

**RKMilGeneral: Generalized Runge-Kutta Milstein Method**

Adaptive step-size Milstein method supporting general (non-diagonal, non-commutative)
noise via iterated stochastic integrals (Levy area approximation).

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive step size
  - **Noise types**: General (scalar, diagonal, non-diagonal, commutative)
  - **SDE interpretation**: Configurable (Ito or Stratonovich)

## Keyword Arguments

  - `interpretation`: `SciMLBase.AlgorithmInterpretation.Ito` (default) or
    `SciMLBase.AlgorithmInterpretation.Stratonovich`
  - `ii_approx`: Iterated integral approximation method (default: `IILevyArea()`)
  - `c`: Truncation parameter for Levy area computation (default: 1)
  - `p`: Truncation level (default: `nothing` for automatic; set `true` with `dt` for manual)
  - `dt`: Time step (used only when `p=true` for manual truncation)

## When to Use

  - For SDEs with general (non-commutative) noise requiring strong order 1.0
  - When adaptive step size is needed with Milstein-type accuracy
  - When Levy area computation is acceptable for accuracy improvement
  - For non-diagonal noise SDEs where RKMil/RKMilCommute are not applicable

## Performance Notes

  - Levy area computation scales with noise dimension
  - Adaptive truncation balances accuracy and efficiency

## References

  - Kloeden, P.E., Platen, E., Numerical Solution of Stochastic Differential Equations. Springer. Berlin Heidelberg (2011)
  - Kastner, F. and Roessler, A., "LevyArea.jl: A Julia package for Levy area computation", arXiv:2201.08424
  - LevyArea.jl: https://github.com/stochastics-uni-luebeck/LevyArea.jl
"""
struct RKMilGeneral{T, TruncationType} <: StochasticDiffEqAdaptiveAlgorithm
    interpretation::SciMLBase.AlgorithmInterpretation.T
    ii_approx::T
    c::Int
    p::TruncationType
end

function RKMilGeneral(;
        interpretation = SciMLBase.AlgorithmInterpretation.Ito,
        ii_approx = IILevyArea(), c = 1, p = nothing, dt = nothing
    )
    gamma = 1 // 1
    p == true && (p = Int(floor(c * dt^(1 // 1 - 2 // 1 * gamma)) + 1))
    return RKMilGeneral{typeof(ii_approx), typeof(p)}(interpretation, ii_approx, c, p)
end

"""
    WangLi3SMil_A()

**WangLi3SMil_A: 3-Stage Milstein Method A (Nonstiff)**

Fixed step-size explicit 3-stage Milstein method with strong and weak order 1.0 for Ito SDEs.

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed step size
  - **Noise types**: Depends on tableau (typically diagonal/scalar)
  - **SDE interpretation**: Ito

## When to Use

  - When fixed step size is preferred
  - For Ito SDEs requiring order 1.0 accuracy
  - Part of WangLi family - compare performance with other variants
  - When computational cost per step is less important than simplicity

## References

  - Wang and Li, "Three-stage stochastic Runge-Kutta methods for stochastic differential equations"
"""
struct WangLi3SMil_A <: StochasticDiffEqAlgorithm end
"""
    WangLi3SMil_B()

**WangLi3SMil_B: 3-Stage Milstein Method B (Nonstiff)**

Alternative 3-stage Milstein method with different stability and accuracy characteristics.

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed step size
  - **Noise types**: Depends on tableau (typically diagonal/scalar)
  - **SDE interpretation**: Ito

## References

  - Wang and Li, "Three-stage stochastic Runge-Kutta methods for stochastic differential equations"
"""
struct WangLi3SMil_B <: StochasticDiffEqAlgorithm end
"""
    WangLi3SMil_C()

**WangLi3SMil_C: 3-Stage Milstein Method C (Nonstiff)**

Third variant in the WangLi 3-stage Milstein family.

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed step size
  - **Noise types**: Depends on tableau (typically diagonal/scalar)
  - **SDE interpretation**: Ito

## References

  - Wang and Li, "Three-stage stochastic Runge-Kutta methods for stochastic differential equations"
"""
struct WangLi3SMil_C <: StochasticDiffEqAlgorithm end
"""
    WangLi3SMil_D()

**WangLi3SMil_D: 3-Stage Milstein Method D (Nonstiff)**

Fourth variant in the WangLi 3-stage Milstein family.

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed step size
  - **Noise types**: Depends on tableau (typically diagonal/scalar)
  - **SDE interpretation**: Ito

## References

  - Wang and Li, "Three-stage stochastic Runge-Kutta methods for stochastic differential equations"
"""
struct WangLi3SMil_D <: StochasticDiffEqAlgorithm end
"""
    WangLi3SMil_E()

**WangLi3SMil_E: 3-Stage Milstein Method E (Nonstiff)**

Fifth variant in the WangLi 3-stage Milstein family.

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed step size
  - **Noise types**: Depends on tableau (typically diagonal/scalar)
  - **SDE interpretation**: Ito

## References

  - Wang and Li, "Three-stage stochastic Runge-Kutta methods for stochastic differential equations"
"""
struct WangLi3SMil_E <: StochasticDiffEqAlgorithm end
"""
    WangLi3SMil_F()

**WangLi3SMil_F: 3-Stage Milstein Method F (Nonstiff)**

Sixth and final variant in the WangLi 3-stage Milstein family.

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed step size
  - **Noise types**: Depends on tableau (typically diagonal/scalar)
  - **SDE interpretation**: Ito

## When to Use (WangLi Family)

  - Compare all variants (A-F) to find best performance for your problem
  - Fixed step applications where step size is predetermined
  - Benchmarking against adaptive methods
  - When Milstein accuracy is needed with explicit fixed steps

## References

  - Wang and Li, "Three-stage stochastic Runge-Kutta methods for stochastic differential equations"
"""
struct WangLi3SMil_F <: StochasticDiffEqAlgorithm end
