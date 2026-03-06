# Abstract types are defined in OrdinaryDiffEqCore and imported via StochasticDiffEq.jl:
# StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
# StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODEAlgorithm,
# StochasticDiffEqRODEAdaptiveAlgorithm, StochasticDiffEqRODECompositeAlgorithm

abstract type StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller} <:
StochasticDiffEqAdaptiveAlgorithm end
abstract type StochasticDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ, Controller} <:
StochasticDiffEqAlgorithm end

abstract type StochasticDiffEqJumpAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpNewtonAdaptiveAlgorithm{
    CS, AD, FDT, ST, CJ, Controller,
} <: StochasticDiffEqJumpAdaptiveAlgorithm end

abstract type StochasticDiffEqJumpDiffusionAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpDiffusionAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm{
    CS, AD, FDT, ST, CJ, Controller,
} <: StochasticDiffEqJumpDiffusionAdaptiveAlgorithm end

abstract type IteratedIntegralApprox end

"""
    IICommutative <: IteratedIntegralApprox

Iterated integral approximation for commutative noise.

This approximation method is used when the noise terms commute, allowing for simplified
computation of iterated stochastic integrals. For commutative noise, only simple stochastic
integrals ∫₀ᵗ dWₛ are needed, and the Lévy area terms vanish.

## When to Use

  - When noise terms satisfy the commutativity condition: `g_i * ∂g_j/∂x_i = g_j * ∂g_i/∂x_j`
  - For diagonal noise SDEs where cross-terms are zero
  - When computational efficiency is more important than handling general non-commutative noise

## Algorithm Properties

  - Assumes zero Lévy area (A_ij = 0 for i ≠ j)
  - Only computes simple stochastic integrals
  - Exact for truly commutative noise structures

!!! warning

    If used on a non-commutative noise problem, this will limit the strong convergence to 0.5,
    regardless of the chosen solver's theoretical order.

## References

  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential Equations", Springer (1992). ISBN: 978-3-540-54062-5. DOI: 10.1007/978-3-662-12616-5
  - Added in PR #459 with the integration of LevyArea.jl package
"""
struct IICommutative <: IteratedIntegralApprox end

"""
    IILevyArea <: IteratedIntegralApprox

Iterated integral approximation using full Lévy area calculation.

This method computes the full Lévy area terms for non-commutative noise, including the
double integrals ∫₀ᵗ ∫₀ˢ dWᵤdWₛ required for higher-order methods. Uses the LevyArea.jl
package which automatically selects optimal algorithms based on problem characteristics.

## When to Use

  - For general non-commutative noise problems
  - When high accuracy is required for methods of strong order > 0.5
  - For problems where noise terms do not commute
  - With RKMilGeneral and other high-order methods

## Algorithm Properties

  - Computes full Lévy area terms A_ij via double stochastic integrals
  - Automatically selects optimal algorithm: Fourier(), Milstein(), Wiktorsson(), or MronRoe()
  - Selection based on Brownian process dimension and step size
  - Required for achieving theoretical convergence rates with non-commutative noise

## Computational Cost

  - More expensive than IICommutative due to Lévy area calculations
  - Cost scales with the number of Brownian motions
  - Optimized algorithm selection minimizes computational overhead

!!! caution

    May introduce bias with adaptive time-stepping methods due to the dependency of
    random number generation on the step size. Use with fixed-step methods when
    maximum accuracy is required.

## References

  - Kastner, F., Rößler, A.: An analysis of approximation algorithms for iterated stochastic
    integrals and a Julia and MATLAB simulation toolbox. Numerical Algorithms 93, 27–66 (2023)
  - Wiktorsson, M. "Joint characteristic function and simultaneous simulation of
    iterated Itô integrals for multiple independent Brownian motions" (2001)
  - Implemented via LevyArea.jl package integration (PR #459)
"""
struct IILevyArea <: IteratedIntegralApprox end

################################################################################

# Basics
"""
EM: Nonstiff Method
The Euler-Maruyama method is the simplest and most fundamental numerical method for solving stochastic differential equations.

## Method Properties

  - **Strong Order**: 0.5 (in the Itô sense)
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed time step only
  - **Noise types**: All forms (diagonal, non-diagonal, scalar, additive, and colored noise)
  - **SDE interpretation**: Itô

## Parameters

  - `split::Bool = true`: Controls step splitting for improved stability with large diffusion eigenvalues

## When to Use

  - First choice for simple SDE problems
  - When computational efficiency is more important than accuracy
  - For problems with all noise types including non-commutative noise
  - When step splitting is needed for stability with large diffusion terms

## Algorithm Description

The method discretizes the SDE:

```math
du = f(u,p,t)dt + g(u,p,t)dW
```

using the scheme:

```math
u_{n+1} = u_n + f(u_n,p,t_n)Δt + g(u_n,p,t_n)ΔW_n
```

When `split=true`, the method uses step splitting which can improve stability for problems with large diffusion eigenvalues.

## References

  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential Equations", Springer (1992). ISBN: 978-3-540-54062-5. DOI: 10.1007/978-3-662-12616-5
"""
struct EM{split} <: StochasticDiffEqAlgorithm end
EM(split = true) = EM{split}()

"""
    SplitEM()

**SplitEM: Split-Step Euler-Maruyama Method (Nonstiff)**

Split-step version of the Euler-Maruyama method that separates the drift and diffusion terms for improved stability.

## Method Properties

  - **Strong Order**: 0.5 (in the Itô sense)
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed time step
  - **Noise types**: All forms (diagonal, non-diagonal, scalar, additive, and colored noise)
  - **SDE interpretation**: Itô

## When to Use

  - When standard EM has stability issues with large diffusion terms
  - Alternative to EM with split=true
  - For problems where operator splitting is natural
  - When drift and diffusion have different timescales

## Algorithm Description

Applies operator splitting to treat drift and diffusion separately:

- Step 1: drift step ``u* = u_n + f(u_n,t_n) Δt``
- Step 2: diffusion step ``u_{n+1} = u* + g(u*,t_n) ΔW_n``

## References

  - Operator splitting methods for SDEs
"""
struct SplitEM <: StochasticDiffEqAlgorithm end
"""
    EulerHeun()

**EulerHeun: Euler-Heun Method (Nonstiff)**

The Euler-Heun method is the Stratonovich analog of the Euler-Maruyama method, providing strong order 0.5 convergence in the Stratonovich sense.

## Method Properties

  - **Strong Order**: 0.5 (in the Stratonovich sense)
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed time step only
  - **Noise types**: All forms (diagonal, non-diagonal, scalar, additive, and colored noise)
  - **SDE interpretation**: Stratonovich

## When to Use

  - When working with Stratonovich SDEs
  - For problems naturally formulated in Stratonovich interpretation
  - When physical interpretation requires Stratonovich calculus
  - As the Stratonovich counterpart to Euler-Maruyama

## Algorithm Description

For Stratonovich SDEs:

```math
du = f(u,p,t)dt + g(u,p,t)∘dW
```

The method uses:

```math
u_{n+1} = u_n + f(u_n,p,t_n)Δt + g(u_n + 0.5*g(u_n,p,t_n)ΔW_n, p, t_n)ΔW_n
```

## Stratonovich vs Itô

  - **EulerHeun**: For Stratonovich SDEs
  - **EM**: For Itô SDEs
  - Conversion between interpretations changes the drift term

## References

  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential Equations", Springer (1992). ISBN: 978-3-540-54062-5. DOI: 10.1007/978-3-662-12616-5
"""
struct EulerHeun <: StochasticDiffEqAlgorithm end
"""
    LambaEM(split=true)

**LambaEM: Adaptive Euler-Maruyama Method (Nonstiff)**

Adaptive time-stepping version of the Euler-Maruyama method with error estimation based on the work of Lamba and Rackauckas.

## Method Properties

  - **Strong Order**: 0.5 (in the Itô sense)
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive with embedded error estimation
  - **Noise types**: All forms (diagonal, non-diagonal, scalar, additive, and colored noise)
  - **SDE interpretation**: Itô

## Parameters

  - `split::Bool = true`: Controls step splitting for improved stability

## When to Use

  - When adaptive time stepping is needed with basic Euler-Maruyama
  - For problems requiring error control without high-order accuracy
  - When computational efficiency and adaptivity are both important
  - For non-commutative noise where higher-order methods aren't applicable

## Algorithm Description

Extends EM with adaptive time stepping using error estimation. The method computes two approximations and uses their difference to estimate local error.

## Error Control

  - Embedded error estimation for adaptive stepping
  - Accepts standard tolerances (abstol, reltol)
  - Automatic step size adjustment

## References

  - Based on error estimation work by Lamba and Rackauckas
"""
struct LambaEM{split} <: StochasticDiffEqAdaptiveAlgorithm end
LambaEM(split = true) = LambaEM{split}()
"""
    LambaEulerHeun()

**LambaEulerHeun: Adaptive Euler-Heun Method (Nonstiff)**

Adaptive time-stepping version of the Euler-Heun method with error estimation for Stratonovich SDEs.

## Method Properties

  - **Strong Order**: 0.5 (in the Stratonovich sense)
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive with embedded error estimation
  - **Noise types**: All forms (diagonal, non-diagonal, scalar, additive, and colored noise)
  - **SDE interpretation**: Stratonovich

## When to Use

  - When adaptive time stepping is needed for Stratonovich SDEs
  - For problems requiring error control in Stratonovich interpretation
  - When computational efficiency and adaptivity are both important
  - For non-commutative noise in Stratonovich formulation

## Algorithm Description

Adaptive version of EulerHeun method with error estimation for automatic step size control.

## Error Control

  - Embedded error estimation for adaptive stepping
  - Standard tolerance control (abstol, reltol)
  - Automatic step size adjustment for Stratonovich problems

## References

  - Error estimation methodology by Lamba, adapted by Rackauckas
"""
struct LambaEulerHeun <: StochasticDiffEqAdaptiveAlgorithm end

"""
    SimplifiedEM()

**SimplifiedEM: Simplified Euler-Maruyama Method (High Weak Order)**

Simplified version of the Euler-Maruyama method optimized for weak convergence with reduced computational cost.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed step size
  - **Noise types**: All forms (diagonal, non-diagonal, scalar, and colored noise)
  - **SDE interpretation**: Itô
  - **Computational cost**: Reduced compared to standard EM

## When to Use

  - Monte Carlo simulations where weak convergence is sufficient
  - Problems where computational efficiency is more important than strong accuracy
  - Large ensemble simulations
  - When only statistical properties (expectations, moments) are needed

## Algorithm Features

  - Simplified implementation reducing computational overhead
  - Maintains weak order 1.0 convergence
  - More efficient than standard EM for weak convergence applications
  - Handles general noise structures

## Weak vs Strong Convergence

  - Optimized for E[f(X_T)] convergence, not pathwise |X_T - X_h| convergence
  - Ideal for computing expectations and statistical properties
  - Less suitable when individual trajectory accuracy is important

## References

  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential Equations", Springer. Berlin Heidelberg (2011)
"""
struct SimplifiedEM <: StochasticDiffEqAlgorithm end

"""
    RKMil(;interpretation=AlgorithmInterpretation.Ito)

**RKMil: Runge-Kutta Milstein Method (Nonstiff)**

Explicit Runge-Kutta discretization of the Milstein method achieving strong order 1.0 convergence for diagonal and scalar noise.

## Method Properties

  - **Strong Order**: 1.0 (for diagonal/scalar noise)
  - **Weak Order**: Depends on tableau
  - **Time stepping**: Adaptive
  - **Noise types**: Diagonal and scalar noise only
  - **SDE interpretation**: Configurable (Itô or Stratonovich)

## Parameters

  - `interpretation`: Choose `AlgorithmInterpretation.Ito` (default) or `AlgorithmInterpretation.Stratonovich`

## When to Use

  - When higher accuracy than Euler methods is needed
  - For diagonal or scalar noise problems
  - When strong order 1.0 convergence is required
  - Alternative to SRI methods for simpler noise structures

## Restrictions

  - **Only works with diagonal or scalar noise**
  - For non-diagonal noise, use RKMilCommute or RKMilGeneral
  - For general noise, use SRI/SRA methods

## Algorithm Description

Implements the Milstein scheme using Runge-Kutta techniques:

```math
du = f(u,t)dt + g(u,t)dW + 0.5 g(u,t) g'(u,t) (dW^2 - dt)
```

where g'(u,t) is the derivative of g with respect to u.

## References

  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential Equations", Springer (1992). ISBN: 978-3-540-54062-5. DOI: 10.1007/978-3-662-12616-5
  - Milstein, G.N., "Numerical Integration of Stochastic Differential Equations"
"""
struct RKMil{interpretation} <: StochasticDiffEqAdaptiveAlgorithm end
RKMil(; interpretation = SciMLBase.AlgorithmInterpretation.Ito) = RKMil{interpretation}()

"""
    RKMilCommute(;interpretation=AlgorithmInterpretation.Ito, ii_approx=IICommutative())

**RKMilCommute: Runge-Kutta Milstein for Commutative Noise (Nonstiff) - Recommended for Commutative Noise**

Explicit Runge-Kutta discretization of the strong order 1.0 Milstein method specialized for **commutative noise** problems.

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: Depends on tableau
  - **Time stepping**: Adaptive (1.5/2.0 error estimate)
  - **Noise types**: Commutative noise (multiple noise sources that commute)
  - **SDE interpretation**: Configurable (Itô or Stratonovich)

## Parameters

  - `interpretation`: Choose `AlgorithmInterpretation.Ito` (default) or `AlgorithmInterpretation.Stratonovich`
  - `ii_approx`: Iterated integral approximation method (default: `IICommutative()`)

## When to Use

  - **Recommended for commutative noise problems**
  - When you have multiple noise sources that satisfy commutativity conditions
  - For multi-dimensional SDEs with commuting noise terms
  - When higher order accuracy than Euler-Maruyama is needed

## Commutative Noise

Applicable when the noise terms satisfy:

```
[g_i, g_j] = g_i(∂g_j/∂x) - g_j(∂g_i/∂x) = 0
```

for all noise terms g_i, g_j.

## Algorithm Description

Extends the Milstein method to handle multiple commutative noise sources efficiently without requiring the full Lévy area computation.

## References

  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential Equations", Springer (1992)
  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential Equations." Springer. Berlin Heidelberg (2011)
"""
struct RKMilCommute{T} <: StochasticDiffEqAdaptiveAlgorithm
    interpretation::SciMLBase.AlgorithmInterpretation.T
    ii_approx::T
end
function RKMilCommute(; interpretation = SciMLBase.AlgorithmInterpretation.Ito, ii_approx = IICommutative())
    return RKMilCommute(interpretation, ii_approx)
end

"""
    RKMilGeneral(;interpretation=AlgorithmInterpretation.Ito, ii_approx=IILevyArea(), c=1, p=nothing, dt=nothing)

**RKMilGeneral: General Runge-Kutta Milstein Method (Nonstiff)**

Explicit Runge-Kutta discretization of the Milstein method for general non-commutative noise problems using Lévy area approximation.

## Method Properties

  - **Strong Order**: 1.0 (for general noise)
  - **Weak Order**: Depends on tableau and Lévy area approximation
  - **Time stepping**: Adaptive
  - **Noise types**: All forms including non-commutative noise
  - **SDE interpretation**: Configurable (Itô or Stratonovich)

## Parameters

  - `interpretation`: Choose `AlgorithmInterpretation.Ito` (default) or `AlgorithmInterpretation.Stratonovich`
  - `ii_approx`: Iterated integral approximation method (default: `IILevyArea()`)
  - `c::Int = 1`: Truncation parameter for Lévy area
  - `p`: Truncation level (computed automatically if `nothing`)
  - `dt`: Used for automatic truncation level computation

## When to Use

  - For general non-commutative noise problems
  - When RKMilCommute is not applicable (noise doesn't commute)
  - When higher accuracy than Euler methods is needed
  - For complex multi-dimensional noise structures

## Lévy Area Approximation

Uses LevyArea.jl for efficient computation of iterated integrals:

  - Automatically selects optimal algorithms
  - Handles truncation for practical computation
  - Supports various approximation strategies

## Computational Cost

  - More expensive than commutative methods
  - Lévy area computation scales with noise dimension
  - Adaptive truncation balances accuracy and efficiency

## References

  - Kloeden, P.E., Platen, E., Numerical Solution of Stochastic Differential Equations. Springer. Berlin Heidelberg (2011)
  - Kastner, F. and Rößler, A., "LevyArea.jl: A Julia package for Lévy area computation", arXiv:2201.08424
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
    γ = 1 // 1
    p == true && (p = Int(floor(c * dt^(1 // 1 - 2 // 1 * γ)) + 1))
    return RKMilGeneral{typeof(ii_approx), typeof(p)}(interpretation, ii_approx, c, p)
end

"""
    WangLi3SMil_A()

**WangLi3SMil_A: 3-Stage Milstein Method A (Nonstiff)**

Fixed step-size explicit 3-stage Milstein method with strong and weak order 1.0 for Itô SDEs.

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed step size
  - **Noise types**: Depends on tableau (typically diagonal/scalar)
  - **SDE interpretation**: Itô

## When to Use

  - When fixed step size is preferred
  - For Itô SDEs requiring order 1.0 accuracy
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
  - **SDE interpretation**: Itô

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
  - **SDE interpretation**: Itô

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
  - **SDE interpretation**: Itô

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
  - **SDE interpretation**: Itô

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
  - **SDE interpretation**: Itô

## When to Use (WangLi Family)

  - Compare all variants (A-F) to find best performance for your problem
  - Fixed step applications where step size is predetermined
  - Benchmarking against adaptive methods
  - When Milstein accuracy is needed with explicit fixed steps

## References

  - Wang and Li, "Three-stage stochastic Runge-Kutta methods for stochastic differential equations"
"""
struct WangLi3SMil_F <: StochasticDiffEqAlgorithm end

#SROCK methods
"""
    SROCK1(;interpretation=AlgorithmInterpretation.Ito, eigen_est=nothing)

**SROCK1: First-Order Stabilized Runge-Kutta Chebyshev Method**

Fixed step size stabilized explicit method designed for mildly stiff SDE problems, particularly effective for parabolic PDEs discretized by method of lines.

## Method Properties

  - **Strong Order**: 0.5 (optimized to 1.0 for scalar/diagonal noise)
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed step size with extended stability
  - **Noise types**: All forms (diagonal, non-diagonal, scalar, additive)
  - **SDE interpretation**: Configurable (Itô or Stratonovich)
  - **Stability**: Extended along negative real axis

## Parameters

  - `interpretation`: Choose `AlgorithmInterpretation.Ito` (default) or `AlgorithmInterpretation.Stratonovich`
  - `eigen_est`: Eigenvalue estimation for stability region (automatic if `nothing`)

## When to Use

  - Parabolic PDEs with stochastic terms
  - Mildly stiff problems where implicit methods are too expensive
  - Large sparse systems from method of lines
  - When stability region extension is more important than high accuracy

## Stability

  - Extends stability region to approximately [-s², 0] where s is number of stages
  - Number of stages chosen based on eigenvalue estimates
  - More efficient than implicit methods for moderate stiffness

## References

  - Chebyshev methods for parabolic stochastic PDEs
  - ROCK methods for stiff problems
"""
struct SROCK1{interpretation, E} <: StochasticDiffEqAlgorithm
    eigen_est::E
end
function SROCK1(; interpretation = SciMLBase.AlgorithmInterpretation.Ito, eigen_est = nothing)
    return SROCK1{interpretation, typeof(eigen_est)}(eigen_est)
end

# Weak Order 2

for Alg in [:SROCK2, :KomBurSROCK2, :SROCKC2]
    @eval begin
        struct $Alg{E} <: StochasticDiffEqAlgorithm
            eigen_est::E
        end
        $Alg(; eigen_est = nothing) = $Alg(eigen_est)
    end
end

@doc """
    SROCK2(;eigen_est=nothing)

**SROCK2: Second-Order Stabilized Runge-Kutta Chebyshev Method**

Second-order stabilized explicit method with weak order 2.0 for mildly stiff SDE problems.

## Method Properties
- **Strong Order**: 1.0
- **Weak Order**: 2.0
- **Time stepping**: Fixed step size with extended stability
- **Stability**: Extended along negative real axis

## When to Use
- When higher accuracy than SROCK1 is needed
- Parabolic PDEs requiring better weak convergence
- Problems where second-order accuracy justifies increased cost

## References
- Second-order ROCK methods for stochastic problems
""" SROCK2

@doc """
    KomBurSROCK2(;eigen_est=nothing)

**KomBurSROCK2: Komori-Burrage Second-Order SROCK Method**

Alternative second-order stabilized method with different coefficients and stability properties.

## Method Properties
- **Strong Order**: 1.0
- **Weak Order**: 2.0
- **Time stepping**: Fixed step size with extended stability
- **Stability**: Extended along negative real axis

## When to Use
- Alternative to SROCK2 with different stability characteristics
- When SROCK2 performance is unsatisfactory
- Benchmarking against other second-order ROCK methods

## References
- Komori and Burrage stabilized methods
""" KomBurSROCK2

@doc """
    SROCKC2(;eigen_est=nothing)

**SROCKC2: Conservative Second-Order SROCK Method**

Conservative second-order stabilized method designed for robust performance.

## Method Properties
- **Strong Order**: 1.0
- **Weak Order**: 2.0
- **Time stepping**: Fixed step size with extended stability
- **Stability**: Conservative stability region, more robust

## When to Use
- When robustness is more important than efficiency
- For difficult problems where other ROCK methods fail
- As a fallback option for problematic cases

## References
- Conservative ROCK methods for stochastic problems
""" SROCKC2

# ROCK stabilization for EM
"""
    SROCKEM(;strong_order_1=true, eigen_est=nothing)

**SROCKEM: ROCK-Stabilized Euler-Maruyama Method**

Fixed step Euler-Maruyama method with first-order ROCK stabilization for handling stiff problems.

## Method Properties

  - **Strong Order**: 1.0 (default) or 0.5 (if `strong_order_1=false`)
  - **Weak Order**: 1.0 (default) or 0.5 (if `strong_order_1=false`)
  - **Time stepping**: Fixed step size with ROCK stabilization
  - **Noise types**: 1-dimensional, diagonal, and multi-dimensional noise
  - **SDE interpretation**: Itô only
  - **Stability**: ROCK stabilization for moderate stiffness

## Parameters

  - `strong_order_1::Bool = true`: Use strong/weak order 1.0 (true) or 0.5 (false)
  - `eigen_est`: Eigenvalue estimation for stability (automatic if `nothing`)

## When to Use

  - Stiff problems where standard EM fails
  - When ROCK stabilization is preferred over full implicit treatment
  - Problems requiring Euler-Maruyama structure with enhanced stability
  - Multi-dimensional stiff SDEs

## Algorithm Description

Combines Euler-Maruyama discretization with ROCK stabilization techniques to extend the stability region without requiring linear solves.

## References

  - ROCK stabilization techniques applied to SDEs
  - Stabilized Euler methods for stiff problems
"""
struct SROCKEM{E} <: StochasticDiffEqAlgorithm
    strong_order_1::Bool
    eigen_est::E
end
SROCKEM(; strong_order_1 = true, eigen_est = nothing) = SROCKEM(strong_order_1, eigen_est)
"""
    SKSROCK(;post_processing=false, eigen_est=nothing)

**SKSROCK: SK-SROCK Stabilized Method**

Fixed step stabilized explicit method for stiff Itô problems with enhanced stability domain and optional post-processing.

## Method Properties

  - **Strong Order**: 0.5 (up to 2.0 with post-processing)
  - **Weak Order**: 1.0 (up to 2.0 with post-processing)
  - **Time stepping**: Fixed step size with enhanced stability
  - **Noise types**: 1-dimensional, diagonal, and multi-dimensional noise
  - **SDE interpretation**: Itô only
  - **Stability**: Better stability domain than SROCK1

## Parameters

  - `post_processing::Bool = false`: Enable post-processing for higher accuracy (experimental)
  - `eigen_est`: Eigenvalue estimation for stability (automatic if `nothing`)

## When to Use

  - Stiff Itô problems requiring better stability than SROCK1
  - Ergodic dynamical systems (with post-processing)
  - Problems where enhanced stability domain is crucial
  - When experimenting with post-processing techniques

## Post-Processing (Experimental)

  - Can achieve order 2 accuracy for ergodic systems
  - Particularly useful for Brownian dynamics
  - Currently under development - use with caution

## Algorithm Features

  - Enhanced stability compared to SROCK1
  - Handles various noise structures
  - Optional post-processing for specialized applications

## References

  - SK-SROCK methods for stochastic problems
  - Post-processing techniques for ergodic systems
"""
struct SKSROCK{E} <: StochasticDiffEqAlgorithm
    post_processing::Bool
    eigen_est::E
end
function SKSROCK(; post_processing = false, eigen_est = nothing)
    return SKSROCK(post_processing, eigen_est)
end
"""
    TangXiaoSROCK2(;version_num=5, eigen_est=nothing)

**TangXiaoSROCK2: Tang-Xiao Second-Order SROCK Method**

Fixed step size stabilized explicit method with multiple variants offering different stability domains.

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: 2.0
  - **Time stepping**: Fixed step size with extended stability
  - **Noise types**: Various (depends on version)
  - **SDE interpretation**: Itô only
  - **Stability**: Version-dependent stability domains

## Parameters

  - `version_num::Int = 5`: Choose version 1-5 with different stability characteristics
  - `eigen_est`: Eigenvalue estimation for stability (automatic if `nothing`)

## When to Use

  - When experimenting with different stability domains
  - Problems requiring weak order 2.0 with fixed steps
  - Benchmarking different ROCK variants
  - **Note**: Currently under development

## Versions

  - Versions 1-5 offer different stability domains
  - Version 5 (default) typically provides good general performance
  - Choose version based on problem-specific stability requirements

## Development Status

  - Method is under active development
  - Use with caution in production code
  - Consider more established ROCK methods for critical applications

## References

  - Tang and Xiao, "Second-order SROCK methods for stochastic problems"
"""
struct TangXiaoSROCK2{E} <: StochasticDiffEqAlgorithm
    version_num::Int
    eigen_est::E
end
function TangXiaoSROCK2(; version_num = 5, eigen_est = nothing)
    return TangXiaoSROCK2(version_num, eigen_est)
end
###############################################################################

# Predictor Corrector
struct PCEuler{T <: Real, F} <: StochasticDiffEqAlgorithm
    theta::T
    eta::T
    ggprime::F
end

"""
    PCEuler(ggprime; theta=1/2, eta=1/2)

**PCEuler: Predictor-Corrector Euler Method (Nonstiff)**

A predictor-corrector variant of the Euler-Maruyama method requiring analytic derivatives
of the diffusion term, with adjustable implicitness parameters for drift-diffusion coupling.

## Method Properties

  - **Strong Order**: 0.5 (in the Itô sense)
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed time step only
  - **Noise types**: General noise with available derivative information
  - **SDE interpretation**: Itô only

## Parameters

  - `ggprime::Function`: The required derivative of the diffusion term

      + For scalar problems: `ggprime = g * ∂g/∂x`
      + For multi-dimensional problems: `ggprime_k = Σ_{j=1...M, i=1...D} g^(j)_i * ∂g^(j)_k/∂x_i`
      + where `g^(j)` corresponds to the noise vector due to the j-th noise channel
      + Must match the in-place/out-of-place specification of the problem

  - `theta::Real = 0.5`: Degree of implicitness in the drift term (default: 0.5)
  - `eta::Real = 0.5`: Degree of implicitness in the diffusion term (default: 0.5)

## When to Use

  - Problems requiring specific drift-diffusion coupling
  - When analytical ggprime function is available
  - Specialized predictor-corrector applications
  - When the derivative `g*g'` provides stability or accuracy benefits

## Algorithm Description

The method uses a predictor-corrector approach with the specific requirement of
computing the derivative of the diffusion coefficient. This additional information
allows for improved handling of drift-diffusion interactions through the adjustable
parameters θ and η.

## Limitations

  - Requires analytical computation of ggprime (cannot be approximated)
  - Fixed time step only (no adaptive versions available)
  - Limited to Itô interpretation

## References

  - Jentzen, A., Kloeden, P.E., "The numerical approximation of stochastic partial differential equations",
    Milan J. Math. 77, 205–244 (2009). https://doi.org/10.1007/s00032-009-0100-0
  - Originally introduced in PR #88 (commit 42e2510) by Tatsuhiro Onodera (2018)

!!! warning

    The derivative `ggprime` must be computed analytically for correctness.# Rossler
    The original paper contains a typo in the definition of ggprime - this
    implementation follows the corrected formulation.
"""
PCEuler(ggprime; theta = 1 / 2, eta = 1 / 2) = PCEuler(theta, eta, ggprime)

################################################################################

# Rossler

"""
    SRA(;tableau=constructSRA1())

**SRA: Configurable Stochastic Runge-Kutta for Additive Noise (Nonstiff)**

Configurable adaptive strong order 1.5 method for additive noise problems with customizable tableaux.

## Method Properties

  - **Strong Order**: 1.5 (for additive noise)
  - **Weak Order**: Depends on tableau (typically 2.0)
  - **Time stepping**: Adaptive
  - **Noise types**: Additive noise (diagonal, non-diagonal, and scalar)
  - **SDE interpretation**: Both Itô and Stratonovich

## Parameters

  - `tableau`: Tableau specification (default: `constructSRA1()`)

## When to Use

  - When custom tableaux are needed for additive noise problems
  - For research and experimentation with SRA methods
  - When default methods don't provide desired characteristics
  - For benchmarking different SRA variants

## Available Tableaux

  - `constructSRA1()`: Default SRA1 tableau
  - Custom tableaux can be constructed for specialized applications

## References

  - Rößler A., "Runge–Kutta Methods for the Strong Approximation of Solutions of Stochastic Differential Equations", SIAM J. Numer. Anal., 48 (3), pp. 922–952. DOI: 10.1137/09076636X.
"""
struct SRA{TabType} <: StochasticDiffEqAdaptiveAlgorithm
    tableau::TabType
end
SRA(; tableau = constructSRA1()) = SRA(tableau)

"""
    SRI(;tableau=constructSRIW1(), error_terms=4)

**SRI: Configurable Stochastic Runge-Kutta for Itô SDEs (Nonstiff)**

Configurable adaptive strong order 1.5 method for diagonal/scalar Itô SDEs with customizable tableaux.

## Method Properties

  - **Strong Order**: 1.5 (for diagonal/scalar noise)
  - **Weak Order**: Depends on tableau (typically 2.0)
  - **Time stepping**: Adaptive
  - **Noise types**: Diagonal and scalar noise only
  - **SDE interpretation**: Itô

## Parameters

  - `tableau`: Tableau specification (default: `constructSRIW1()`)
  - `error_terms::Int = 4`: Number of error terms for adaptive stepping

## When to Use

  - When custom tableaux are needed for diagonal/scalar problems
  - For research and experimentation with SRI methods
  - When default methods don't provide desired characteristics
  - For benchmarking different SRI variants

## Available Tableaux

  - `constructSRIW1()`: Default SRIW1 tableau
  - Custom tableaux can be constructed for specialized applications

## References

  - Rößler A., "Runge–Kutta Methods for the Strong Approximation of Solutions of Stochastic Differential Equations", SIAM J. Numer. Anal., 48 (3), pp. 922–952. DOI: 10.1137/09076636X.
"""
struct SRI{TabType} <: StochasticDiffEqAdaptiveAlgorithm
    tableau::TabType
    error_terms::Int
end
SRI(; tableau = constructSRIW1(), error_terms = 4) = SRI(tableau, error_terms)

"""
    SRIW1()

**SRIW1: Stochastic Runge-Kutta W1 Method (Nonstiff)**

Adaptive stochastic Runge-Kutta method with strong order 1.5 and weak order 2.0 for diagonal/scalar Itô SDEs.

## Method Properties

  - **Strong Order**: 1.5 (for diagonal/scalar noise)
  - **Weak Order**: 2.0
  - **Time stepping**: Adaptive
  - **Noise types**: Diagonal and scalar noise only
  - **SDE interpretation**: Itô

## When to Use

  - Standard choice for diagonal/scalar Itô SDEs
  - When proven theoretical properties are important
  - Alternative to SOSRI when stability optimization is not needed
  - For problems requiring exactly weak order 2.0

## Algorithm Features

  - Embedded error estimation for adaptive stepping
  - Well-established theoretical foundation
  - Good balance of accuracy and efficiency

## References

  - Rößler A., "Runge–Kutta Methods for the Strong Approximation of Solutions of Stochastic Differential Equations", SIAM J. Numer. Anal., 48 (3), pp. 922–952. DOI: 10.1137/09076636X.
"""
struct SRIW1 <: StochasticDiffEqAdaptiveAlgorithm end

"""
    SRIW2()

**SRIW2: Stochastic Runge-Kutta W2 Method (Nonstiff)**

Adaptive stochastic Runge-Kutta method with strong order 1.5 and weak order 3.0 for diagonal/scalar Itô SDEs.

## Method Properties

  - **Strong Order**: 1.5 (for diagonal/scalar noise)
  - **Weak Order**: 3.0
  - **Time stepping**: Adaptive
  - **Noise types**: Diagonal and scalar noise only
  - **SDE interpretation**: Itô

## When to Use

  - When weak order 3.0 convergence is required
  - For Monte Carlo simulations needing high weak accuracy
  - Problems where weak convergence is more important than strong
  - When computational cost per step is acceptable for higher weak order

## Algorithm Features

  - Highest weak order in the SRI family
  - More expensive per step than SRIW1
  - Excellent for statistical calculations and expectations

## References

  - Rößler A., "Runge–Kutta Methods for the Strong Approximation of Solutions of Stochastic Differential Equations", SIAM J. Numer. Anal., 48 (3), pp. 922–952. DOI: 10.1137/09076636X.
"""
struct SRIW2 <: StochasticDiffEqAdaptiveAlgorithm end
"""
    SOSRI()

**SOSRI: Stability-Optimized SRI Method (Nonstiff) - Recommended**

The Stability-Optimized Stochastic Runge-Kutta method. This is the **recommended method** for general-purpose solving of diagonal/scalar Itô SDEs.

## Method Properties

  - **Strong Order**: 1.5 (for diagonal/scalar noise)
  - **Weak Order**: 2.0
  - **Time stepping**: Adaptive
  - **Noise types**: Diagonal and scalar noise only
  - **SDE interpretation**: Itô
  - **Stability**: Optimized for high tolerances and robust to mild stiffness

## When to Use

  - **Recommended as first choice** for diagonal/scalar Itô SDEs
  - When high accuracy is required (strong order 1.5)
  - For problems with mild stiffness
  - When using high tolerances (method is stable)
  - For most general SDE applications

## Algorithm Description

SOSRI is a stability-optimized version of the SRI methods with specially chosen coefficients to improve stability properties. It provides excellent performance for the most common class of SDE problems.

## Restrictions

  - Only works with diagonal or scalar noise
  - For non-diagonal noise, use other methods like `RKMilCommute` or `LambaEM`

## References

  - Rößler A., "Runge–Kutta Methods for the Strong Approximation of Solutions of Stochastic Differential Equations", SIAM J. Numer. Anal., 48 (3), pp. 922–952
"""
struct SOSRI <: StochasticDiffEqAdaptiveAlgorithm end
"""
    SOSRI2()

**SOSRI2: Alternative Stability-Optimized SRI Method (Nonstiff)**

Alternative stability-optimized adaptive strong order 1.5 method with different stability characteristics than SOSRI.

## Method Properties

  - **Strong Order**: 1.5 (for diagonal/scalar noise)
  - **Weak Order**: 2.0
  - **Time stepping**: Adaptive
  - **Noise types**: Diagonal and scalar noise only
  - **SDE interpretation**: Itô
  - **Stability**: Optimized for high tolerances and robust to stiffness

## When to Use

  - Alternative to SOSRI with different stability properties
  - When SOSRI performance is unsatisfactory
  - For benchmarking stability-optimized methods
  - Problems requiring different stability characteristics

## Algorithm Features

  - Different stability optimization than SOSRI
  - May perform better on certain problem types
  - Maintains high tolerance robustness

## References

  - Stability-optimized SRI methods
"""
struct SOSRI2 <: StochasticDiffEqAdaptiveAlgorithm end

"""
    SRA1()

**SRA1: Stochastic Runge-Kutta A1 Method (Nonstiff)**

Adaptive strong order 1.5 method for additive Itô and Stratonovich SDEs with weak order 2.

## Method Properties

  - **Strong Order**: 1.5 (for additive noise)
  - **Weak Order**: 2.0
  - **Time stepping**: Adaptive
  - **Noise types**: Additive noise (diagonal, non-diagonal, and scalar)
  - **SDE interpretation**: Both Itô and Stratonovich

## When to Use

  - Standard choice for additive noise problems
  - When proven theoretical properties are important
  - Alternative to SOSRA when stability optimization is not needed
  - For both Itô and Stratonovich problems with additive noise

## Additive Noise Structure

Specialized for SDEs of the form:

```math
du = f(u,p,t) dt + σ(p,t) dW
```

where diffusion σ doesn't depend on solution u.

## References

  - Rößler A., "Runge–Kutta Methods for the Strong Approximation of Solutions of Stochastic Differential Equations", SIAM J. Numer. Anal., 48 (3), pp. 922–952. DOI:10.1137/09076636X
"""
struct SRA1 <: StochasticDiffEqAdaptiveAlgorithm end

"""
    SRA2()

**SRA2: Stochastic Runge-Kutta A2 Method (Nonstiff)**

Alternative adaptive strong order 1.5 method for additive noise problems with different coefficients.

## Method Properties

  - **Strong Order**: 1.5 (for additive noise)
  - **Weak Order**: 2.0
  - **Time stepping**: Adaptive
  - **Noise types**: Additive noise (diagonal, non-diagonal, and scalar)
  - **SDE interpretation**: Both Itô and Stratonovich

## When to Use

  - Alternative to SRA1 with different stability/accuracy characteristics
  - When SRA1 performance is unsatisfactory
  - For benchmarking different SRA variants
  - Research and comparison studies

## References

  - Rößler A., "Runge–Kutta Methods for the Strong Approximation of Solutions of Stochastic Differential Equations", SIAM J. Numer. Anal., 48 (3), pp. 922–952. DOI:10.1137/09076636X
"""
struct SRA2 <: StochasticDiffEqAdaptiveAlgorithm end

"""
    SRA3()

**SRA3: Stochastic Runge-Kutta A3 Method (Nonstiff)**

Adaptive strong order 1.5 method for additive noise problems with weak order 3.

## Method Properties

  - **Strong Order**: 1.5 (for additive noise)
  - **Weak Order**: 3.0
  - **Time stepping**: Adaptive
  - **Noise types**: Additive noise (non-diagonal and scalar)
  - **SDE interpretation**: Both Itô and Stratonovich

## When to Use

  - When weak order 3.0 convergence is required for additive noise
  - For Monte Carlo simulations needing highest weak accuracy
  - Problems where weak convergence dominates computational cost
  - When computational cost per step is acceptable for higher weak order

## Restrictions

  - **Does not handle diagonal additive noise** (use SRA1/SRA2 instead)
  - Limited to non-diagonal and scalar additive noise structures

## Algorithm Features

  - Highest weak order in the SRA family
  - More expensive per step than SRA1/SRA2
  - Excellent for statistical calculations requiring high weak accuracy

## References

  - Rößler A., "Runge–Kutta Methods for the Strong Approximation of Solutions of Stochastic Differential Equations", SIAM J. Numer. Anal., 48 (3), pp. 922–952. DOI:10.1137/09076636X
"""
struct SRA3 <: StochasticDiffEqAdaptiveAlgorithm end
"""
    SOSRA()

**SOSRA: Stability-Optimized SRA Method (Nonstiff) - Optimal for Additive Noise**

Stability-optimized adaptive Stochastic Runge-Kutta method for additive noise problems. This is the **optimal choice for additive noise SDEs**.

## Method Properties

  - **Strong Order**: 1.5 (for additive noise)
  - **Weak Order**: 2.0
  - **Time stepping**: Adaptive
  - **Noise types**: Additive noise (diagonal, non-diagonal, and scalar)
  - **SDE interpretation**: Both Itô and Stratonovich
  - **Stability**: Optimized for high tolerances and robust to stiffness

## When to Use

  - **Optimal choice** for additive noise problems: du = f(u,p,t)dt + σ dW
  - When the diffusion term is independent of the solution u
  - For problems requiring high accuracy with additive noise
  - When using high tolerances (method is stable)
  - For both Itô and Stratonovich interpretations

## Algorithm Description

SOSRA is a stability-optimized version of the SRA (Stochastic Runge-Kutta for Additive noise) methods. It exploits the special structure of additive noise to achieve better performance and stability.

## Additive Noise Structure

Specialized for SDEs of the form:

```math
du = f(u,p,t) dt + σ(t) dW
```

where the diffusion σ does not depend on the solution u.

## References

  - Rößler A., "Runge–Kutta Methods for the Strong Approximation of Solutions of Stochastic Differential Equations", SIAM J. Numer. Anal., 48 (3), pp. 922–952
"""
struct SOSRA <: StochasticDiffEqAdaptiveAlgorithm end
"""
    SOSRA2()

**SOSRA2: Stability-Optimized SRA Method Version 2 (Nonstiff)**

Alternative stability-optimized adaptive SRA method for additive noise problems.

## Method Properties

  - **Strong Order**: 1.5 (for additive noise)
  - **Weak Order**: 2.0
  - **Time stepping**: Adaptive
  - **Noise types**: Additive noise (diagonal, non-diagonal, and scalar)
  - **SDE interpretation**: Both Itô and Stratonovich
  - **Stability**: Optimized for high tolerances and robust to stiffness

## When to Use

  - Alternative to SOSRA for additive noise problems
  - Different stability characteristics may be preferred for specific problems
  - When SOSRA performance is unsatisfactory

## References

  - Rößler A., "Runge–Kutta Methods for the Strong Approximation of Solutions of Stochastic Differential Equations", SIAM J. Numer. Anal., 48 (3), pp. 922–952
"""
struct SOSRA2 <: StochasticDiffEqAdaptiveAlgorithm end

################################################################################

# Rossler second order for weak approx.

"""
    DRI1()

**DRI1: Debrabant-Rößler Implicit Method (High Weak Order)**

Adaptive high-order method optimized for weak convergence with minimized error constants. Excellent for Monte Carlo simulations and moment calculations.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0 (optimized with minimized error constants)
  - **Deterministic Order**: 3.0 (when noise = 0)
  - **Time stepping**: Adaptive
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive noise)
  - **SDE interpretation**: Itô

## When to Use

  - **Optimal for weak convergence** requirements
  - Monte Carlo simulations where statistical properties matter most
  - Computing expectations, moments, and probability distributions
  - When weak accuracy is more important than pathwise accuracy
  - For problems requiring diverse noise types

## Weak vs Strong Convergence

  - **Weak convergence**: Convergence of expectations E[f(X_T)]
  - **Strong convergence**: Pathwise convergence |X_T - X_T^h|
  - DRI1 prioritizes weak convergence with optimized error constants

## Algorithm Features

  - Minimized error constants for better practical performance
  - Handles complex noise structures including non-commuting terms
  - Adaptive time stepping for efficiency

## References

  - Debrabant, K. and Rößler A., "Families of efficient second order Runge–Kutta methods for the weak approximation of Itô stochastic differential equations", Applied Numerical Mathematics 59, pp. 582–594 (2009). DOI: 10.1016/j.apnum.2008.03.012
"""
struct DRI1 <: StochasticDiffEqAdaptiveAlgorithm end

"""
    DRI1NM()

**DRI1NM: Debrabant-Rößler Implicit Non-Mixing Method (High Weak Order)**

Specialized version of DRI1 for non-mixing diagonal and scalar additive noise problems.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0 (optimized with minimized error constants)
  - **Deterministic Order**: 3.0 (when noise = 0)
  - **Time stepping**: Adaptive
  - **Noise types**: Non-mixing diagonal and scalar additive noise
  - **SDE interpretation**: Itô

## When to Use

  - Non-mixing diagonal problems: `du[k] = f(u[k]) dt + σ[k] dW[k]`
  - Scalar additive noise problems
  - When DRI1 is too general/expensive for the problem structure
  - Monte Carlo simulations with special structure

## Non-Mixing Diagonal Structure

Optimized for problems where:

```
du[1] = f₁(u[1])dt + σ₁ dW[1]
du[2] = f₂(u[2])dt + σ₂ dW[2]
...
```

Each component depends only on itself (no coupling).

## Algorithm Advantages

  - More efficient than general DRI1 for structured problems
  - Exploits special structure for better performance
  - Maintains weak order 2.0 with minimized constants

## References

  - Debrabant, K. and Rößler A., "Families of efficient second order Runge–Kutta methods for the weak approximation of Itô stochastic differential equations", Applied Numerical Mathematics 59, pp. 582–594 (2009). DOI: 10.1016/j.apnum.2008.03.012.
"""
struct DRI1NM <: StochasticDiffEqAdaptiveAlgorithm end

"""
    RI1()

**RI1: Rößler Implicit Method 1 (High Weak Order)**

Adaptive weak order 2.0 method for Itô SDEs with deterministic order 3.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 3.0 (when noise = 0)
  - **Time stepping**: Adaptive
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive)
  - **SDE interpretation**: Itô

## When to Use

  - General weak convergence problems
  - Monte Carlo simulations with various noise structures
  - When weak order 2.0 is sufficient
  - Alternative to DRI1 with different characteristics

## References

  - Rößler A., "Second Order Runge–Kutta Methods for Itô Stochastic Differential Equations", SIAM J. Numer. Anal., 47, pp. 1713-1738 (2009). DOI: 10.1137/060673308.
"""
struct RI1 <: StochasticDiffEqAdaptiveAlgorithm end

"""
    RI3()

**RI3: Rößler Implicit Method 3 (High Weak Order)**

Alternative adaptive weak order 2.0 method with different stability characteristics.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 3.0 (when noise = 0)
  - **Time stepping**: Adaptive
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive)
  - **SDE interpretation**: Itô

## When to Use

  - Alternative to RI1 with different characteristics
  - When RI1 performance is unsatisfactory
  - Benchmarking different weak order 2.0 methods

## References

  - Rößler A., "Second Order Runge–Kutta Methods for Itô Stochastic Differential Equations", SIAM J. Numer. Anal., 47, pp. 1713-1738 (2009). DOI: 10.1137/060673308
"""
struct RI3 <: StochasticDiffEqAdaptiveAlgorithm end

"""
    RI5()

**RI5: Rößler Implicit Method 5 (High Weak Order)**

Another variant in the RI family of weak order 2.0 methods.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 3.0 (when noise = 0)
  - **Time stepping**: Adaptive
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive)
  - **SDE interpretation**: Itô

## When to Use

  - Part of RI family comparison studies
  - When other RI methods don't provide desired characteristics
  - Research applications requiring different RI variants

## References

  - Rößler A., "Second Order Runge–Kutta Methods for Itô Stochastic Differential Equations", SIAM J. Numer. Anal., 47, pp. 1713-1738 (2009). DOI: 10.1137/060673308
"""
struct RI5 <: StochasticDiffEqAdaptiveAlgorithm end

"""
    RI6()

**RI6: Rößler Implicit Method 6 (High Weak Order)**

Final method in the RI family with deterministic order 2 (lower than other RI methods).

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 2.0 (when noise = 0)
  - **Time stepping**: Adaptive
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive)
  - **SDE interpretation**: Itô

## When to Use

  - When lower deterministic order is acceptable
  - Potentially more efficient than RI1/RI3/RI5
  - Completing RI family comparisons

## Algorithm Features

  - Lower deterministic order may reduce computational cost
  - Still maintains weak order 2.0 for stochastic problems
  - Final variant in the comprehensive RI family

## References

  - Rößler A., "Second Order Runge–Kutta Methods for Itô Stochastic Differential Equations", SIAM J. Numer. Anal., 47, pp. 1713-1738 (2009). DOI: 10.1137/060673308
"""
struct RI6 <: StochasticDiffEqAdaptiveAlgorithm end

"""
    RDI1WM()

**RDI1WM: Runge-Kutta Debrabant Implicit 1 Weak Method (High Weak Order)**

Fixed step method with weak order 1.0 for Itô SDEs.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 1.0
  - **Deterministic Order**: 2.0 (when noise = 0)
  - **Time stepping**: Fixed step size
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive)
  - **SDE interpretation**: Itô

## When to Use

  - Fixed step applications where step size is predetermined
  - When weak order 1.0 is sufficient
  - Simpler alternative to higher-order weak methods
  - Baseline for comparing higher-order methods

## References

  - Debrabant, K. and Rößler A., "Classification of Stochastic Runge–Kutta Methods for the Weak Approximation of Stochastic Differential Equations", Mathematics and Computers in Simulation 77, pp. 408-420 (2008). DOI: 10.1016/j.matcom.2007.04.016
"""
struct RDI1WM <: StochasticDiffEqAlgorithm end

"""
    RDI2WM()

**RDI2WM: Runge-Kutta Debrabant Implicit 2 Weak Method (High Weak Order)**

Adaptive weak order 2.0 method for Itô SDEs with deterministic order 2.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 2.0 (when noise = 0)
  - **Time stepping**: Adaptive
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive)
  - **SDE interpretation**: Itô

## When to Use

  - Weak order 2.0 problems with adaptive stepping
  - Alternative to DRI1 and RI methods
  - When deterministic order 2.0 is sufficient
  - Monte Carlo simulations requiring adaptive control

## References

  - Debrabant, K. and Rößler A., "Classification of Stochastic Runge–Kutta Methods for the Weak Approximation of Stochastic Differential Equations", Mathematics and Computers in Simulation 77, pp. 408-420 (2008). DOI: 10.1016/j.matcom.2007.04.016.
"""
struct RDI2WM <: StochasticDiffEqAdaptiveAlgorithm end

"""
    RDI3WM()

**RDI3WM: Runge-Kutta Debrabant Implicit 3 Weak Method (High Weak Order)**

Adaptive weak order 2.0 method with higher deterministic order 3.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 3.0 (when noise = 0)
  - **Time stepping**: Adaptive
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive)
  - **SDE interpretation**: Itô

## When to Use

  - When both weak order 2.0 and deterministic order 3.0 are needed
  - Problems with significant deterministic components
  - Alternative to DRI1 with different characteristics
  - High accuracy requirements for both stochastic and deterministic parts

## References

  - Debrabant, K. and Rößler A., "Classification of Stochastic Runge–Kutta Methods for the Weak Approximation of Stochastic Differential Equations", Mathematics and Computers in Simulation 77, pp. 408-420 (2008). DOI: 10.1016/j.matcom.2007.04.016.
"""
struct RDI3WM <: StochasticDiffEqAdaptiveAlgorithm end

"""
    RDI4WM()

**RDI4WM: Runge-Kutta Debrabant Implicit 4 Weak Method (High Weak Order)**

Fourth variant in the RDI family with weak order 2.0 and deterministic order 3.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 3.0 (when noise = 0)
  - **Time stepping**: Adaptive
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive)
  - **SDE interpretation**: Itô

## When to Use

  - Final alternative in the RDI family
  - When other RDI methods don't provide desired performance
  - Completing comprehensive RDI method comparisons
  - Research applications requiring all RDI variants

## References

  - Debrabant, K. and Rößler A., "Classification of Stochastic Runge–Kutta Methods for the Weak Approximation of Stochastic Differential Equations", Mathematics and Computers in Simulation 77, pp. 408-420 (2008). DOI: 10.1016/j.matcom.2007.04.016.
"""
struct RDI4WM <: StochasticDiffEqAdaptiveAlgorithm end

"""
    W2Ito1()

**W2Ito1: Wang-Tang-Xiao Weak Order 2 Method (High Weak Order)**

Efficient weak second-order method for Itô SDEs with adaptive stepping.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 3.0 (when noise = 0)
  - **Time stepping**: Adaptive
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive)
  - **SDE interpretation**: Itô

## When to Use

  - Modern efficient weak order 2.0 method
  - When computational efficiency is important for weak convergence
  - Alternative to older weak order 2.0 methods
  - Monte Carlo simulations requiring good performance

## Algorithm Features

  - Designed for computational efficiency
  - Good balance of accuracy and cost for weak problems
  - More recent development than classical methods

## References

  - Tang, X., & Xiao, A., "Efficient weak second-order stochastic Runge–Kutta methods for Itô stochastic differential equations", BIT Numerical Mathematics, 57, 241-260 (2017). DOI: 10.1007/s10543-016-0618-9
"""
struct W2Ito1 <: StochasticDiffEqAdaptiveAlgorithm end

# Stratonovich sense

"""
    RS1()

**RS1: Rößler Stratonovich Method 1 (High Weak Order)**

Fixed step weak order 2.0 method specifically designed for Stratonovich SDEs.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 2.0 (when noise = 0)
  - **Time stepping**: Fixed step size
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive)
  - **SDE interpretation**: Stratonovich

## When to Use

  - Stratonovich SDEs requiring weak order 2.0
  - Fixed step applications with predetermined step size
  - Problems naturally formulated in Stratonovich interpretation
  - When physical interpretation requires Stratonovich calculus

## Stratonovich Interpretation

Optimized for SDEs in Stratonovich form:

```math
du = f(u,t)dt + g(u,t)∘dW
```

where ∘ denotes Stratonovich integration.

## References

  - Rößler A., "Second order Runge–Kutta methods for Stratonovich stochastic differential equations", BIT Numerical Mathematics 47, pp. 657-680 (2007) DOI: 10.1007/s10543-007-0130-3.
"""
struct RS1 <: StochasticDiffEqAlgorithm end

"""
    RS2()

**RS2: Rößler Stratonovich Method 2 (High Weak Order)**

Alternative fixed step weak order 2.0 method for Stratonovich SDEs with higher deterministic order.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 3.0 (when noise = 0)
  - **Time stepping**: Fixed step size
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive)
  - **SDE interpretation**: Stratonovich

## When to Use

  - Stratonovich SDEs with significant deterministic components
  - When higher deterministic accuracy than RS1 is needed
  - Fixed step applications requiring better deterministic performance
  - Benchmarking against RS1

## Algorithm Features

  - Higher deterministic order than RS1
  - May be more expensive per step than RS1
  - Better for problems with large deterministic components

## References

  - Rößler A., "Second order Runge–Kutta methods for Stratonovich stochastic differential equations", BIT Numerical Mathematics 47, pp. 657-680 (2007). DOI: 10.1007/s10543-007-0130-3
"""
struct RS2 <: StochasticDiffEqAlgorithm end

"""
    PL1WM()

**PL1WM: Platen Weak Method 1 (High Weak Order)**

Fixed step weak order 2.0 method from the classical Kloeden-Platen textbook.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 2.0 (when noise = 0)
  - **Time stepping**: Fixed step size
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive)
  - **SDE interpretation**: Itô

## When to Use

  - Classical reference implementation for weak order 2.0
  - Fixed step applications with predetermined step size
  - Educational purposes and textbook examples
  - Baseline comparison for more advanced methods

## Algorithm Features

  - Well-established classical method
  - Simple implementation
  - Standard reference from foundational SDE literature

## References

  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential Equations", Springer. Berlin Heidelberg (2011). ISBN 978-3-540-54062-5. DOI: 10.1007/978-3-662-12616-5.
"""
struct PL1WM <: StochasticDiffEqAlgorithm end

"""
    PL1WMA()

**PL1WMA: Platen Weak Method 1 Additive (High Weak Order)**

Specialized version of PL1WM optimized for additive noise problems.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 2.0 (when noise = 0)
  - **Time stepping**: Fixed step size
  - **Noise types**: Additive noise only
  - **SDE interpretation**: Itô

## When to Use

  - Additive noise problems with fixed step size
  - When PL1WM is too general for additive structure
  - Classical reference for additive noise weak methods
  - Educational and benchmarking purposes

## Additive Noise Structure

Specialized for SDEs of the form:

```math
du = f(u,t)dt + σ(t) dW
```

where diffusion σ doesn't depend on solution u.

## Algorithm Features

  - More efficient than PL1WM for additive problems
  - Classical foundation method
  - Simplified implementation for additive case

## References

  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential Equations", Springer. Berlin Heidelberg (2011). ISBN 978-3-540-54062-5. DOI: 10.1007/978-3-662-12616-5.
"""
struct PL1WMA <: StochasticDiffEqAlgorithm end

"""
    NON()

NON: High Weak Order Method
Fixed step weak order 2.0 for Stratonovich SDEs (deterministic order 4).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.

## References

  - Komori, Y., Weak second-order stochastic Runge–Kutta methods for non-commutative stochastic differential equations, Journal of Computational and Applied Mathematics 206, pp. 158 – 173 (2007). DOI: 10.1016/j.cam.2006.06.006.
"""
struct NON <: StochasticDiffEqAlgorithm end

"""
    COM()

**COM: Commutative Stratonovich Method (High Weak Order)**

Fixed step method optimized for commutative Stratonovich SDEs.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: Depends on implementation
  - **Time stepping**: Fixed step size
  - **Noise types**: Commutative noise only
  - **SDE interpretation**: Stratonovich

## When to Use

  - Commutative Stratonovich SDEs
  - When noise terms satisfy commutativity conditions
  - More efficient alternative to NON for commutative cases
  - Fixed step applications with commutative structure

## Commutative Noise

Optimized for Stratonovich SDEs where:

```
[g_i, g_j] = g_i(∂g_j/∂x) - g_j(∂g_i/∂x) = 0
```

for all noise terms.

## Algorithm Features

  - More efficient than NON for commutative cases
  - Exploits commutativity for computational savings
  - Specialized for Stratonovich interpretation

## References

  - Komori, Y., "Weak order stochastic Runge–Kutta methods for commutative stochastic differential equations", Journal of Computational and Applied Mathematics 203, pp. 57 – 79 (2007). Bibcode: 2007JCoAM.203...57K
"""
struct COM <: StochasticDiffEqAlgorithm end

"""
    NON2()

**NON2: Enhanced Non-commutative Stratonovich Method (High Weak Order)**

Improved version of the NON method with enhanced efficiency for non-commutative Stratonovich SDEs.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Time stepping**: Fixed step size
  - **Noise types**: Non-commutative noise
  - **SDE interpretation**: Stratonovich

## When to Use

  - Enhanced version of NON with better efficiency
  - Non-commutative Stratonovich SDEs requiring improved performance
  - When NON is too expensive or inefficient
  - Modern alternative to classical NON method

## Algorithm Features

  - More efficient than original NON method
  - Maintains weak order 2.0 convergence
  - Enhanced computational techniques

## References

  - Komori, Y., & Burrage, K., "Supplement: Efficient weak second order stochastic Runge–Kutta methods for non-commutative Stratonovich stochastic differential equations", Journal of computational and applied mathematics, 235(17), pp. 5326-5329 (2011)
"""
struct NON2 <: StochasticDiffEqAlgorithm end

"""
    SIEA()

**SIEA: Stochastic Improved Euler A Method (High Weak Order)**

Stochastic generalization of the improved Euler method for Itô SDEs.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 2.0 (when noise = 0)
  - **Time stepping**: Fixed step size
  - **Noise types**: Diagonal and scalar additive noise
  - **SDE interpretation**: Itô

## When to Use

  - Fixed step applications with diagonal/scalar additive noise
  - When stochastic version of improved Euler is desired
  - Educational purposes (connection to classical methods)
  - Baseline for Tocino-Vigo-Aguiar method comparisons

## Algorithm Features

  - Based on classical improved Euler method
  - Specialized for additive noise structures
  - Simple and well-understood foundation

## References

  - Tocino, A. and Vigo-Aguiar, J., "Weak Second Order Conditions for Stochastic Runge-Kutta Methods", SIAM Journal on Scientific Computing 24, pp. 507-523 (2002). DOI: 10.1137/S1064827501387814.
"""
struct SIEA <: StochasticDiffEqAlgorithm end

"""
    SMEA()

**SMEA: Stochastic Modified Euler A Method (High Weak Order)**

Stochastic generalization of the modified Euler method for Itô SDEs.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 2.0 (when noise = 0)
  - **Time stepping**: Fixed step size
  - **Noise types**: Diagonal and scalar additive noise
  - **SDE interpretation**: Itô

## When to Use

  - Fixed step applications with diagonal/scalar additive noise
  - When stochastic version of modified Euler is desired
  - Alternative to SIEA with different characteristics
  - Educational and comparison purposes

## Algorithm Features

  - Based on classical modified Euler method
  - Different approach than SIEA for same problem class
  - Specialized for additive noise structures

## References

  - Tocino, A. and Vigo-Aguiar, J., "Weak Second Order Conditions for Stochastic Runge-Kutta Methods", SIAM Journal on Scientific Computing 24, pp. 507-523 (2002). DOI: 10.1137/S1064827501387814.
"""
struct SMEA <: StochasticDiffEqAlgorithm end

"""
    SIEB()

**SIEB: Stochastic Improved Euler B Method (High Weak Order)**

Alternative stochastic generalization of the improved Euler method.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 2.0 (when noise = 0)
  - **Time stepping**: Fixed step size
  - **Noise types**: Diagonal and scalar additive noise
  - **SDE interpretation**: Itô

## When to Use

  - Alternative to SIEA with different coefficients
  - Fixed step applications requiring different stability properties
  - Comparing different improved Euler generalizations
  - When SIEA performance is unsatisfactory

## Algorithm Features

  - Variant B of stochastic improved Euler approach
  - Different coefficients than SIEA
  - May have different stability or accuracy characteristics

## References

  - Tocino, A. and Vigo-Aguiar, J., "Weak Second Order Conditions for Stochastic Runge-Kutta Methods", SIAM Journal on Scientific Computing 24, pp. 507-523 (2002). DOI: 10.1137/S1064827501387814.
"""
struct SIEB <: StochasticDiffEqAlgorithm end

"""
    SMEB()

**SMEB: Stochastic Modified Euler B Method (High Weak Order)**

Alternative stochastic generalization of the modified Euler method.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 2.0 (when noise = 0)
  - **Time stepping**: Fixed step size
  - **Noise types**: Diagonal and scalar additive noise
  - **SDE interpretation**: Itô

## When to Use

  - Alternative to SMEA with different coefficients
  - Fixed step applications requiring different characteristics
  - Completing Tocino-Vigo-Aguiar method family comparisons
  - When SMEA performance is unsatisfactory

## Algorithm Features

  - Variant B of stochastic modified Euler approach
  - Different coefficients than SMEA
  - Completes the family of Tocino-Vigo-Aguiar methods

## References

  - Tocino, A. and Vigo-Aguiar, J., "Weak Second Order Conditions for Stochastic Runge-Kutta Methods", SIAM Journal on Scientific Computing 24, pp. 507-523 (2002). DOI: 10.1137/S1064827501387814
"""
struct SMEB <: StochasticDiffEqAlgorithm end

################################################################################

# IIF

"""
    IIF1M(;nlsolve=NLSOLVEJL_SETUP())

**IIF1M: Integrating Factor Method 1 (Semi-Linear)**

First-order integrating factor method for semi-linear SDEs with stiff linear parts.

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed or adaptive
  - **Problem type**: Semi-linear SDEs with stiff linear components
  - **Treatment**: Exponential integrator approach

## Parameters

  - `nlsolve`: Nonlinear solver configuration

## When to Use

  - Semi-linear SDEs: du = (L*u + N(u))dt + g(u)dW where L is stiff linear operator
  - Problems amenable to integrating factor techniques
  - When exponential integrators are appropriate
  - Stiff linear parts with nonlinear perturbations

## Algorithm Description

Applies integrating factor exp(L*t) to handle stiff linear part exactly while treating nonlinear parts numerically.

## References

  - Integrating factor methods for stiff SDEs
"""
struct IIF1M{F} <: StochasticDiffEqAlgorithm
    nlsolve::F
end
IIF1M(; nlsolve = NLSOLVEJL_SETUP()) = IIF1M{typeof(nlsolve)}(nlsolve)

"""
    IIF2M(;nlsolve=NLSOLVEJL_SETUP())

**IIF2M: Integrating Factor Method 2 (Semi-Linear)**

Second-order integrating factor method for semi-linear SDEs.

## Method Properties

  - **Strong Order**: 2.0
  - **Weak Order**: 2.0
  - **Time stepping**: Fixed or adaptive
  - **Problem type**: Semi-linear SDEs with stiff linear components
  - **Treatment**: Higher-order exponential integrator

## Parameters

  - `nlsolve`: Nonlinear solver configuration

## When to Use

  - When higher accuracy than IIF1M is needed
  - Semi-linear problems requiring second-order convergence
  - More expensive but more accurate than IIF1M

## References

  - Higher-order integrating factor methods for SDEs
"""
struct IIF2M{F} <: StochasticDiffEqAlgorithm
    nlsolve::F
end
IIF2M(; nlsolve = NLSOLVEJL_SETUP()) = IIF2M{typeof(nlsolve)}(nlsolve)

"""
    IIF1Mil(;nlsolve=NLSOLVEJL_SETUP())

**IIF1Mil: Integrating Factor Milstein Method (Semi-Linear)**

Integrating factor method combined with Milstein correction for semi-linear SDEs.

## Method Properties

  - **Strong Order**: 1.0 (with Milstein correction)
  - **Weak Order**: 1.0
  - **Time stepping**: Fixed or adaptive
  - **Problem type**: Semi-linear SDEs with stiff linear components
  - **Treatment**: Exponential integrator with Milstein correction

## Parameters

  - `nlsolve`: Nonlinear solver configuration

## When to Use

  - Semi-linear SDEs requiring Milstein-type accuracy
  - When both stiff linear treatment and higher-order stochastic accuracy are needed
  - Alternative to IIF1M with enhanced stochastic treatment

## References

  - Integrating factor methods with Milstein correction
"""
struct IIF1Mil{F} <: StochasticDiffEqAlgorithm
    nlsolve::F
end
IIF1Mil(; nlsolve = NLSOLVEJL_SETUP()) = IIF1Mil{typeof(nlsolve)}(nlsolve)

################################################################################

# SDIRK
"""
    ImplicitEM(;chunk_size=0, autodiff=true, diff_type=Val{:central},
               standardtag=Val{true}(), concrete_jac=nothing, precs=DEFAULT_PRECS,
               linsolve=nothing, nlsolve=NLNewton(), extrapolant=:constant,
               theta=1, symplectic=false, new_jac_conv_bound=1e-3, 
               controller=:Predictive)

**ImplicitEM: Implicit Euler-Maruyama Method (Stiff)**

Drift-implicit version of the Euler-Maruyama method with theta-method treatment of the drift term.

## Method Properties

  - **Strong Order**: 0.5 (Itô sense)
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive (1.0/1.5 heuristic)
  - **Noise types**: All forms (non-diagonal, scalar, colored noise)
  - **SDE interpretation**: Itô
  - **Implicit treatment**: Drift term only (diffusion remains explicit)

## Parameters

  - `theta::Real = 1`: Implicitness parameter (0=explicit, 1=fully implicit, 0.5=trapezoidal)
  - `symplectic::Bool = false`: When `true` and `theta=0.5`, uses symplectic implicit midpoint
  - Linear/nonlinear solver options via `linsolve` and `nlsolve`

## When to Use

  - For mildly stiff SDEs where drift term causes stability issues
  - When explicit methods require very small time steps
  - As a robust fallback for difficult problems
  - When all noise types need to be supported

## Theta Method Variants

  - `theta = 0`: Explicit Euler (not recommended, use `EM` instead)
  - `theta = 0.5`: Trapezoidal rule (second order accurate for deterministic part)
  - `theta = 1`: Backward Euler (default, maximum stability)

## Symplectic Option

When `symplectic=true` and `theta=0.5`, the method preserves the symplectic structure in distribution for appropriate problems.

## Algorithm Description

Treats the SDE ``du = f(u,t)dt + g(u,t)dW`` using:

```math
u_{n+1} = u_n + θ f(u_{n+1},t_{n+1}) dt + (1-θ) f(u_n,t_n) dt + g(u_n,t_n) dW_n
```

## References

  - Standard implicit methods adapted for SDEs
"""
struct ImplicitEM{CS, AD, F, F2, P, FDT, ST, CJ, T2, Controller} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller}
    linsolve::F
    nlsolve::F2
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
end
function ImplicitEM(;
        chunk_size = 0, autodiff = true, diff_type = Val{:central},
        standardtag = Val{true}(), concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    return ImplicitEM{
        chunk_size, autodiff,
        typeof(linsolve), typeof(nlsolve), typeof(precs), diff_type,
        SciMLBase._unwrap_val(standardtag),
        SciMLBase._unwrap_val(concrete_jac),
        typeof(new_jac_conv_bound), controller,
    }(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant, new_jac_conv_bound, symplectic
    )
end

STrapezoid(; kwargs...) = ImplicitEM(; theta = 1 / 2, kwargs...)
SImplicitMidpoint(; kwargs...) = ImplicitEM(; theta = 1 / 2, symplectic = true, kwargs...)
"""
    ImplicitEulerHeun(;chunk_size=0, autodiff=true, diff_type=Val{:central},
                      standardtag=Val{true}(), concrete_jac=nothing, precs=DEFAULT_PRECS,
                      linsolve=nothing, nlsolve=NLNewton(), extrapolant=:constant,
                      theta=1, symplectic=false, new_jac_conv_bound=1e-3, 
                      controller=:Predictive)

**ImplicitEulerHeun: Implicit Euler-Heun Method (Stiff)**

Drift-implicit version of the Euler-Heun method for Stratonovich SDEs with stiff drift terms.

## Method Properties

  - **Strong Order**: 0.5 (Stratonovich sense)
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive (1.0/1.5 heuristic)
  - **Noise types**: All forms (non-diagonal, scalar, colored noise)
  - **SDE interpretation**: Stratonovich
  - **Implicit treatment**: Drift term only (diffusion remains explicit)

## Parameters

  - `theta::Real = 1`: Implicitness parameter (0=explicit, 1=fully implicit, 0.5=trapezoidal)
  - `symplectic::Bool = false`: When `true` and `theta=1`, uses symplectic implicit midpoint
  - Linear/nonlinear solver options via `linsolve` and `nlsolve`

## When to Use

  - Stiff Stratonovich SDEs where drift term causes stability issues
  - When working in Stratonovich interpretation with stiff dynamics
  - Alternative to ImplicitEM for Stratonovich problems
  - When all noise types need to be supported in Stratonovich form

## Theta Method Variants

  - `theta = 0.5`: Trapezoidal rule (default, good accuracy/stability balance)
  - `theta = 1`: Backward Euler (maximum stability)

## Symplectic Option

When `symplectic=true` and `theta=1`, preserves symplectic structure in distribution for appropriate Stratonovich problems.

## References

  - Implicit methods for stiff SDEs in Stratonovich interpretation
"""
struct ImplicitEulerHeun{CS, AD, F, P, FDT, ST, CJ, N, T2, Controller} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller}
    linsolve::F
    nlsolve::N
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
end
function ImplicitEulerHeun(;
        chunk_size = 0, autodiff = true, diff_type = Val{:central},
        standardtag = Val{true}(), concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    return ImplicitEulerHeun{
        chunk_size, autodiff,
        typeof(linsolve), typeof(precs), diff_type,
        SciMLBase._unwrap_val(standardtag),
        SciMLBase._unwrap_val(concrete_jac),
        typeof(nlsolve),
        typeof(new_jac_conv_bound), controller,
    }(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic
    )
end

"""
    ImplicitRKMil(;chunk_size=0, autodiff=true, diff_type=Val{:central},
                  standardtag=Val{true}(), concrete_jac=nothing, precs=DEFAULT_PRECS,
                  linsolve=nothing, nlsolve=NLNewton(), extrapolant=:constant,
                  theta=1, symplectic=false, new_jac_conv_bound=1e-3, 
                  controller=:Predictive, interpretation=AlgorithmInterpretation.Ito)

**ImplicitRKMil: Implicit Runge-Kutta Milstein Method (Stiff)**

Drift-implicit Runge-Kutta Milstein method achieving order 1.0 for stiff problems with diagonal/scalar noise.

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: Depends on tableau
  - **Time stepping**: Adaptive (1.5/2.0 heuristic)
  - **Noise types**: Diagonal and scalar noise only
  - **SDE interpretation**: Configurable (Itô or Stratonovich)
  - **Implicit treatment**: Drift term only (diffusion remains explicit)

## Parameters

  - `theta::Real = 1`: Implicitness parameter (0.5=trapezoidal, 1=backward Euler)
  - `symplectic::Bool = false`: When `true` and `theta=0.5`, uses symplectic implicit midpoint
  - `interpretation`: Choose `AlgorithmInterpretation.Ito` (default) or `AlgorithmInterpretation.Stratonovich`
  - Linear/nonlinear solver options via `linsolve` and `nlsolve`

## When to Use

  - Stiff problems requiring higher accuracy than ImplicitEM
  - When strong order 1.0 is needed with implicit stability
  - Diagonal or scalar noise problems with stiff drift
  - Alternative to SKenCarp for non-additive noise

## Restrictions

  - **Only works with diagonal or scalar noise**
  - For non-diagonal noise, use ISSEM/ISSEulerHeun
  - For additive noise, prefer SKenCarp

## Algorithm Features

  - Higher order accuracy than ImplicitEM
  - Milstein correction for improved strong convergence
  - Configurable interpretation (Itô/Stratonovich)

## References

  - Implicit Milstein methods for stiff SDEs
"""
struct ImplicitRKMil{CS, AD, F, P, FDT, ST, CJ, N, T2, Controller, interpretation} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller}
    linsolve::F
    nlsolve::N
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
end
function ImplicitRKMil(;
        chunk_size = 0, autodiff = true, diff_type = Val{:central},
        standardtag = Val{true}(), concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive, interpretation = SciMLBase.AlgorithmInterpretation.Ito
    )
    return ImplicitRKMil{
        chunk_size, autodiff,
        typeof(linsolve), typeof(precs), diff_type,
        SciMLBase._unwrap_val(standardtag),
        SciMLBase._unwrap_val(concrete_jac),
        typeof(nlsolve), typeof(new_jac_conv_bound),
        controller, interpretation,
    }(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic
    )
end
"""
    ISSEM(;chunk_size=0, autodiff=true, diff_type=Val{:central},
          standardtag=Val{true}(), concrete_jac=nothing, precs=DEFAULT_PRECS,
          linsolve=nothing, nlsolve=NLNewton(), extrapolant=:constant,
          theta=1, symplectic=false, new_jac_conv_bound=1e-3, 
          controller=:Predictive)

**ISSEM: Implicit Split-Step Euler-Maruyama Method (Stiff)**

Fully implicit split-step method for handling stiffness in both drift and diffusion terms.

## Method Properties

  - **Strong Order**: 0.5 (Itô sense)
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive (1.0/1.5 heuristic)
  - **Noise types**: All forms (non-diagonal, scalar, colored noise)
  - **SDE interpretation**: Itô
  - **Implicit treatment**: Both drift and diffusion terms (fully implicit)

## Parameters

  - `theta::Real = 1`: Implicitness parameter for drift term
  - `symplectic::Bool = false`: When `true` and `theta=0.5`, uses symplectic implicit midpoint
  - Linear/nonlinear solver options via `linsolve` and `nlsolve`

## When to Use

  - **Recommended for stiff Itô problems with large noise terms**
  - When both drift and diffusion cause stability issues
  - Problems where ImplicitEM and ImplicitRKMil are insufficient
  - Fully stiff SDEs requiring implicit treatment of everything

## Algorithm Description

Applies implicit treatment to both drift and diffusion using split-step approach:

- Step 1: Handle drift implicitly
- Step 2: Handle diffusion implicitly

## Fully Implicit Features

  - Can handle stiffness in both drift and diffusion
  - More expensive than drift-only implicit methods
  - Most robust for extremely stiff problems
  - Requires solving nonlinear systems for both terms

## References

  - Split-step implicit methods for fully stiff SDEs
"""
struct ISSEM{CS, AD, F, P, FDT, ST, CJ, N, T2, Controller} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller}
    linsolve::F
    nlsolve::N
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
end
function ISSEM(;
        chunk_size = 0, autodiff = true, diff_type = Val{:central},
        standardtag = Val{true}(), concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    return ISSEM{
        chunk_size, autodiff,
        typeof(linsolve), typeof(precs), diff_type,
        SciMLBase._unwrap_val(standardtag),
        SciMLBase._unwrap_val(concrete_jac),
        typeof(nlsolve),
        typeof(new_jac_conv_bound), controller,
    }(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic
    )
end
"""
    ISSEulerHeun(;chunk_size=0, autodiff=true, diff_type=Val{:central},
                 standardtag=Val{true}(), concrete_jac=nothing, precs=DEFAULT_PRECS,
                 linsolve=nothing, nlsolve=NLNewton(), extrapolant=:constant,
                 theta=1, symplectic=false, new_jac_conv_bound=1e-3, 
                 controller=:Predictive)

**ISSEulerHeun: Implicit Split-Step Euler-Heun Method (Stiff)**

Fully implicit split-step method for Stratonovich SDEs with stiffness in both drift and diffusion terms.

## Method Properties

  - **Strong Order**: 0.5 (Stratonovich sense)
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive (1.0/1.5 heuristic)
  - **Noise types**: All forms (non-diagonal, scalar, colored noise)
  - **SDE interpretation**: Stratonovich
  - **Implicit treatment**: Both drift and diffusion terms (fully implicit)

## Parameters

  - `theta::Real = 1`: Implicitness parameter for drift term
  - `symplectic::Bool = false`: When `true` and `theta=0.5`, uses symplectic implicit midpoint
  - Linear/nonlinear solver options via `linsolve` and `nlsolve`

## When to Use

  - **Recommended for stiff Stratonovich problems with large noise terms**
  - When both drift and diffusion cause stability issues in Stratonovich form
  - Stratonovich problems where ImplicitEulerHeun is insufficient
  - Fully stiff Stratonovich SDEs

## Algorithm Description

Stratonovich analog of ISSEM with fully implicit treatment of both drift and diffusion terms using split-step approach.

## Fully Implicit Features

  - Handles stiffness in both drift and diffusion for Stratonovich SDEs
  - Most expensive but most robust for Stratonovich stiff problems
  - Requires solving nonlinear systems for both drift and diffusion

## References

  - Split-step implicit methods for fully stiff Stratonovich SDEs
"""
struct ISSEulerHeun{CS, AD, F, P, FDT, ST, CJ, N, T2, Controller} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller}
    linsolve::F
    nlsolve::N
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
end
function ISSEulerHeun(;
        chunk_size = 0, autodiff = true, diff_type = Val{:central},
        standardtag = Val{true}(), concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    return ISSEulerHeun{
        chunk_size, autodiff,
        typeof(linsolve), typeof(precs), diff_type,
        SciMLBase._unwrap_val(standardtag),
        SciMLBase._unwrap_val(concrete_jac),
        typeof(nlsolve), typeof(new_jac_conv_bound), controller,
    }(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic
    )
end
"""
    SKenCarp(;chunk_size=0, autodiff=true, diff_type=Val{:central}, 
             standardtag=Val{true}(), concrete_jac=nothing, precs=DEFAULT_PRECS,
             linsolve=nothing, nlsolve=NLNewton(), smooth_est=true, 
             extrapolant=:min_correct, new_jac_conv_bound=1e-3, 
             controller=:Predictive, ode_error_est=true)

**SKenCarp: Stochastic KenCarp Method (Stiff) - Highly Recommended for Stiff Problems**

Adaptive L-stable drift-implicit method with strong order 1.5. **Highly recommended** for stiff problems with additive noise.

## Method Properties

  - **Strong Order**: 1.5 (for additive noise)
  - **Weak Order**: 2.0
  - **Time stepping**: Adaptive
  - **Noise types**: Additive noise (diagonal, non-diagonal, and scalar)
  - **SDE interpretation**: Both Itô and Stratonovich
  - **Stability**: L-stable (excellent for stiff problems)
  - **Implicit**: Drift-implicit (handles stiffness in drift term)

## When to Use

  - **Highly recommended** for stiff additive noise problems
  - When the drift term f(u,p,t) is stiff
  - For problems requiring high accuracy with stiff dynamics
  - When implicit treatment of the drift is necessary for stability
  - Best choice for stiff problems with additive noise structure

## Algorithm Description

SKenCarp applies implicit treatment to the drift term while keeping the diffusion explicit. This provides excellent stability for stiff SDEs with additive noise.

## Stiffness and Stability

  - L-stable: Excellent for stiff problems
  - Handles large negative eigenvalues in the drift term
  - Maintains accuracy while providing stability

## Configuration Options

  - Linear solver options via `linsolve` parameter
  - Nonlinear solver options via `nlsolve` parameter
  - Jacobian computation control via `autodiff` and related parameters
  - Step size control via `controller` parameter

## References

  - Based on KenCarp methods from OrdinaryDiffEq.jl
  - Adapted for stochastic problems with additive noise
"""
struct SKenCarp{CS, AD, F, P, FDT, ST, CJ, N, T2, Controller} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller}
    linsolve::F
    nlsolve::N
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    new_jac_conv_bound::T2
    ode_error_est::Bool
end

function SKenCarp(;
        chunk_size = 0, autodiff = true, diff_type = Val{:central},
        standardtag = Val{true}(), concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :min_correct,
        new_jac_conv_bound = 1.0e-3, controller = :Predictive,
        ode_error_est = true
    )
    return SKenCarp{
        chunk_size, autodiff, typeof(linsolve), typeof(precs), diff_type,
        SciMLBase._unwrap_val(standardtag), SciMLBase._unwrap_val(concrete_jac),
        typeof(nlsolve), typeof(new_jac_conv_bound), controller,
    }(
        linsolve, nlsolve, precs, smooth_est, extrapolant, new_jac_conv_bound,
        ode_error_est
    )
end

################################################################################

# Implicit Weak Order 2 Methods

"""
    IRI1(;chunk_size=0, autodiff=true, diff_type=Val{:central},
         standardtag=Val{true}(), concrete_jac=nothing, precs=DEFAULT_PRECS,
         linsolve=nothing, nlsolve=NLNewton(), extrapolant=:constant,
         theta=1, new_jac_conv_bound=1e-3, controller=:Predictive)

**IRI1: Implicit Rößler 1 Method (Stiff, High Weak Order)**

Drift-implicit version of the RI1 weak order 2.0 method for stiff Itô SDEs with multiplicative noise.

## Method Properties

  - **Strong Order**: Not optimized for strong convergence
  - **Weak Order**: 2.0
  - **Deterministic Order**: 3.0 (when noise = 0)
  - **Time stepping**: Adaptive
  - **Noise types**: All forms (diagonal, non-diagonal, non-commuting, scalar additive)
  - **SDE interpretation**: Itô
  - **Implicit treatment**: Drift term only (diffusion remains explicit)

## Parameters

  - `theta::Real = 1`: Implicitness parameter (0=explicit, 1=fully implicit, 0.5=trapezoidal)
  - Linear/nonlinear solver options via `linsolve` and `nlsolve`

## When to Use

  - **Stiff SDEs with multiplicative noise requiring weak order 2.0**
  - When SKenCarp cannot be used (non-additive noise)
  - Monte Carlo simulations where drift causes stability issues
  - When weak convergence is sufficient but drift stability is needed
  - Alternative to explicit RI1 when drift is stiff

## Algorithm Description

IRI1 applies the theta-method implicitization to the drift stages of the RI1 weak order 2
stochastic Runge-Kutta method. The diffusion terms remain explicit, making this method
suitable for problems where the drift is stiff but the diffusion does not cause stability issues.

## Theta Method Variants

  - `theta = 0.5`: Trapezoidal rule (good accuracy/stability balance)
  - `theta = 1`: Backward Euler (maximum stability, default)

## Comparison with Other Methods

  - **vs RI1**: IRI1 adds implicit drift treatment for stiff problems
  - **vs SKenCarp**: IRI1 handles multiplicative noise, not just additive
  - **vs ImplicitEM**: IRI1 achieves weak order 2.0 instead of 1.0
  - **vs ISSEM**: IRI1 has higher weak order but only drift-implicit

## References

  - Rößler A., "Second Order Runge–Kutta Methods for Itô Stochastic Differential Equations", SIAM J. Numer. Anal., 47, pp. 1713-1738 (2009). DOI: 10.1137/060673308
  - Debrabant, K., Rößler, A., "Diagonally drift-implicit Runge-Kutta methods of weak order one and two for Itô SDEs and stability analysis", Applied Numerical Mathematics 59(3–4), 595–607 (2009). DOI: 10.1016/j.apnum.2008.03.011
  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential Equations", Springer (1992). Chapter 15: Explicit and Implicit Weak Approximations.
"""
struct IRI1{CS, AD, F, F2, P, FDT, ST, CJ, T2, Controller} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller}
    linsolve::F
    nlsolve::F2
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
end
function IRI1(;
        chunk_size = 0, autodiff = true, diff_type = Val{:central},
        standardtag = Val{true}(), concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    return IRI1{
        chunk_size, autodiff,
        typeof(linsolve), typeof(nlsolve), typeof(precs), diff_type,
        SciMLBase._unwrap_val(standardtag),
        SciMLBase._unwrap_val(concrete_jac),
        typeof(new_jac_conv_bound), controller,
    }(
        linsolve, nlsolve, precs,
        theta,
        extrapolant, new_jac_conv_bound
    )
end

################################################################################

# Jumps

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
struct TauLeaping <: StochasticDiffEqJumpAdaptiveAlgorithm end

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

  - Cao, Y., Gillespie, D.T., Petzold, L.R., "Efficient step size selection for the tau-leaping method"# Etc.
"""
struct CaoTauLeaping <: StochasticDiffEqJumpAdaptiveAlgorithm end

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

################################################################################

# Etc.

"""
    StochasticCompositeAlgorithm(algs, choice_function)

**StochasticCompositeAlgorithm: Multi-Method Composite Algorithm**

Composite algorithm that automatically switches between multiple SDE solvers based on problem characteristics.

## Method Properties

  - **Approach**: Multi-method solving with automatic switching
  - **Adaptivity**: Changes methods during integration
  - **Flexibility**: Combines strengths of different algorithms

## Parameters

  - `algs::Tuple`: Tuple of algorithms to switch between
  - `choice_function::Function`: Function determining which algorithm to use

## When to Use

  - Problems with changing characteristics during integration
  - When different regions require different solution approaches
  - Combining methods for different regimes (e.g., stiff/nonstiff)
  - When no single method is optimal for entire domain

## Choice Function

The choice_function(integrator) should return an integer indicating which algorithm from algs to use:

```julia
function choice_function(integrator)
    if stiff_region(integrator.u, integrator.t)
        return 1  # Use first algorithm (e.g., implicit)
    else
        return 2  # Use second algorithm (e.g., explicit)
    end
end
```

## Algorithm Features

  - Automatic method switching during integration
  - Maintains continuity across method transitions
  - Combines computational efficiency with robustness
  - Can handle complex multi-scale problems

## References

  - Composite algorithm methodology for SDEs
"""
struct StochasticCompositeAlgorithm{T, F} <: StochasticDiffEqCompositeAlgorithm
    algs::T
    choice_function::F
end

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

const SplitSDEAlgorithms = Union{IIF1M, IIF2M, IIF1Mil, SKenCarp, SplitEM}

"""
    BAOAB(;gamma=1.0, scale_noise=true)

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
