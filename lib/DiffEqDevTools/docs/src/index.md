# DiffEqDevTools

Public API for convergence testing, work-precision benchmarks, and Runge–Kutta
tableau utilities.

## Convergence

```@docs
DiffEqDevTools.ConvergenceSimulation
DiffEqDevTools.test_convergence
DiffEqDevTools.analyticless_test_convergence
DiffEqDevTools.appxtrue
DiffEqDevTools.TestSolution
```

## Benchmarks

```@docs
DiffEqDevTools.Shootout
DiffEqDevTools.ShootoutSet
DiffEqDevTools.WorkPrecision
DiffEqDevTools.WorkPrecisionSet
DiffEqDevTools.get_sample_errors
```

## Tableaus

```@docs
DiffEqDevTools.stability_region
DiffEqDevTools.imaginary_stability_interval
DiffEqDevTools.check_tableau
DiffEqDevTools.deduce_Butcher_tableau
DiffEqDevTools.residual_order_condition
```
