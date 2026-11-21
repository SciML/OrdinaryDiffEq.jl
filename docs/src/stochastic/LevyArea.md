# Lévy Area and Iterated Stochastic Integrals

`StochasticDiffEqLevyArea` provides Fourier-series approximations of Lévy areas
and iterated stochastic integrals. It supports direct random sampling,
pre-generated coefficients for deterministic reuse, and Brownian-path
reconstruction within a time step.

## Algorithms

```@docs
StochasticDiffEqLevyArea.AbstractIteratedIntegralAlgorithm
StochasticDiffEqLevyArea.Fourier
StochasticDiffEqLevyArea.Milstein
StochasticDiffEqLevyArea.Wiktorsson
StochasticDiffEqLevyArea.MronRoe
```

## Error norms

```@docs
StochasticDiffEqLevyArea.AbstractErrorNorm
StochasticDiffEqLevyArea.MaxL2
StochasticDiffEqLevyArea.FrobeniusL2
```

## Coefficients and iterated integrals

```@docs
StochasticDiffEqLevyArea.LevyAreaCoefficients
StochasticDiffEqLevyArea.coefficient_length
StochasticDiffEqLevyArea.generate_coefficients
StochasticDiffEqLevyArea.levyarea
StochasticDiffEqLevyArea.terms_needed
StochasticDiffEqLevyArea.optimal_algorithm
StochasticDiffEqLevyArea.iterated_integrals
StochasticDiffEqLevyArea.reconstruct_path
StochasticDiffEqLevyArea.iterated_integrals_subinterval
```
