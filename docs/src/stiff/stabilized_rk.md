# Stabilized Runge-Kutta Methods (Runge-Kutta-Chebyshev)

## Explicit Stabilized Runge-Kutta Methods

Explicit stabilized methods utilize an upper bound on the spectral radius of the Jacobian.
Users can supply an upper bound by specifying the keyword argument `eigen_est`, for example

```julia
`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`
```

The methods `ROCK2` and `ROCK4` also include keyword arguments `min_stages` and `max_stages`,
which specify upper and lower bounds on the adaptively chosen number of stages for stability.

```@docs
ROCK2
ROCK4
SERK2
ESERK4
ESERK5
RKC
```

## Implicit Stabilized Runge-Kutta Methods

```@docs
IRKC
```
