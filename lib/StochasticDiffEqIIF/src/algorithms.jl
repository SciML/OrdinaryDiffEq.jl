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
