# StochasticDiffEqLevyArea.jl

*Iterated Stochastic Integrals with Separated Coefficient Generation*

This package is a successor to [LevyArea.jl](https://github.com/stochastics-uni-luebeck/LevyArea.jl) by Kastner and Rößler, carrying over all algorithms and extending them with support for pre-generated Fourier coefficients, explicit RNG control, Brownian path reconstruction, and sub-interval iterated integral computation. These additions enable consistent iterated integrals across adaptive step-size changes and convergence testing in SDE solvers.

**Users of this package should cite the original LevyArea.jl work:**

> Kastner, F. and Rößler, A., "LevyArea.jl: A Julia package for Lévy area computation", arXiv: [2201.08424](https://arxiv.org/abs/2201.08424)

The underlying algorithms are described in:

> Milstein, G. N., "Numerical integration of stochastic differential equations." Vol. 313. Springer Science & Business Media, 1995. DOI: [10.1007/978-94-015-8455-5](https://doi.org/10.1007/978-94-015-8455-5)

> Wiktorsson, M., "Joint Characteristic Function and Simultaneous Simulation of Iterated Itô Integrals for Multiple Independent Brownian Motions." The Annals of Applied Probability 11.2 (2001), pp. 470-487. DOI: [10.1214/aoap/1015345301](https://doi.org/10.1214/aoap/1015345301)

> Mrongowius, J. and Rößler, A., "On the Approximation and Simulation of Iterated Stochastic Integrals and the Corresponding Lévy Areas in Terms of a Multidimensional Brownian Motion." Stochastic Analysis and Applications (2021), pp. 1-29. arXiv: [2101.09542](https://arxiv.org/abs/2101.09542). DOI: [10.1080/07362994.2021.1922291](https://doi.org/10.1080/07362994.2021.1922291)

## What's New vs LevyArea.jl

All four algorithms from LevyArea.jl (Fourier, Milstein, Wiktorsson, MronRoe) are carried over with bit-identical output when using the default RNG. The new APIs are:

- **`LevyAreaCoefficients`**: Store pre-generated Fourier coefficients separately from the area computation, enabling deterministic and reproducible Levy area given the same coefficients.
- **`generate_coefficients(m, n, alg, rng)`**: Generate coefficients from an explicit RNG, isolating randomness from the global RNG state.
- **`levyarea(W, n, alg, coeffs)`**: Compute Levy area deterministically from pre-generated coefficients.
- **`levyarea(W, n, alg; rng=...)`**: All algorithms accept an explicit `rng` keyword argument.
- **`reconstruct_path(dW, h, coeffs, t_points)`**: Reconstruct the Brownian motion at arbitrary time points within [0, h] from the Karhunen-Loeve expansion using stored Fourier coefficients.
- **`iterated_integrals_subinterval(dW, h, coeffs, t_start, t_end)`**: Compute iterated integrals over sub-intervals by reconstructing the path and using Riemann sums.

## Usage

The LevyArea.jl API is fully preserved:

```julia
using StochasticDiffEqLevyArea

m = 5
h = 0.01
W = sqrt(h) * randn(m)

# Default usage (same as LevyArea.jl)
II = iterated_integrals(W, h)

# With explicit algorithm
II = iterated_integrals(W, h; alg=MronRoe())

# With explicit precision
II = iterated_integrals(W, h, 0.05)
```

New coefficient-based workflow:

```julia
using Random

# Generate and store Fourier coefficients
rng = Xoshiro(42)
n = terms_needed(m, h, h^(3/2), MronRoe(), MaxL2())
coeffs = generate_coefficients(m, n, MronRoe(), rng)

# Deterministic computation from stored coefficients
II = levyarea(W / sqrt(h), n, MronRoe(), coeffs)

# Reconstruct the Brownian path at arbitrary times
t_points = [0.0, h/4, h/2, 3h/4, h]
W_values = reconstruct_path(W, h, coeffs, t_points)

# Compute iterated integrals over a sub-interval
II_sub = iterated_integrals_subinterval(W, h, coeffs, 0.0, h/2)
```

## Available Algorithms

| Algorithm | Convergence Order | Description |
|-----------|-------------------|-------------|
| `Fourier()` | O(n^{-1/2}) | Basic Fourier expansion |
| `Milstein()` | O(n^{-1/2}) | Fourier with simple rest approximation |
| `Wiktorsson()` | O(n^{-1}) | Fourier with tail sum approximation |
| `MronRoe()` | O(n^{-1}) | Improved tail sum (recommended) |
