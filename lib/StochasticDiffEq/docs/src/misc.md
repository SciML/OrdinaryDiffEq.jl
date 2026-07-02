# Miscellaneous Methods

This page covers specialized methods for particular types of problems or applications.

## Composite Algorithms

### StochasticCompositeAlgorithm - Multi-Method Solving

```@docs
StochasticCompositeAlgorithm
```

## RODE Methods (Random ODEs)

### RandomEM - Random Euler Method

```@docs
RandomEM
```

### RandomHeun - Random Heun Method

```@docs
RandomHeun
```

### RandomTamedEM - Tamed Random Euler

```@docs
RandomTamedEM
```

## Langevin Dynamics

### BAOAB - Langevin Integrator

```@docs
BAOAB
```

## Predictor-Corrector Methods

### PCEuler - Predictor-Corrector Euler

```@docs
PCEuler
```

## Integro-Integral-Form (IIF) Methods

### IIF1M, IIF2M, IIF1Mil - IIF Methods

```@docs
IIF1M
```

```@docs
IIF2M
```

```@docs
IIF1Mil
```

## Simplified Methods

### SimplifiedEM - Simplified Euler-Maruyama

```@docs
SimplifiedEM
```

## When to Use Miscellaneous Methods

### StochasticCompositeAlgorithm:

  - When problem characteristics change during integration
  - Combining methods for different regimes
  - Automatic method switching based on conditions

### RODE Methods:

  - Random ordinary differential equations
  - Problems with random parameters but no Brownian motion
  - Uncertainty quantification applications

### BAOAB:

  - Molecular dynamics simulations
  - Langevin equations with specific structure
  - When preserving equilibrium distributions is important

### IIF Methods:

  - Semi-linear problems with stiff linear parts
  - Problems amenable to integrating factor techniques
  - When exponential integrators are appropriate

### PCEuler:

  - Problems requiring specific drift-diffusion coupling
  - When analytical ggprime function is available
  - Specialized predictor-corrector applications

These methods serve specific niches in stochastic computation and may be optimal for particular problem structures.
