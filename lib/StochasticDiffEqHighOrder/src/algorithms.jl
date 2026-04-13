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
