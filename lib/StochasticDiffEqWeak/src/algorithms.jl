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
