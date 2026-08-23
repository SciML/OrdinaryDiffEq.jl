# DiffEqDevTools

Public API for convergence testing, work-precision benchmarks, and Runge–Kutta
tableau utilities.

## Convergence

```@docs
DiffEqDevTools.ConvergenceSimulation
DiffEqDevTools.ConvergenceTrajectory
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

## Tagging and comparison helpers

```@docs
DiffEqDevTools.get_tags
DiffEqDevTools.unique_tags
DiffEqDevTools.filter_by_tags
DiffEqDevTools.exclude_by_tags
DiffEqDevTools.merge_wp_sets
DiffEqDevTools.available_errors
DiffEqDevTools.wp_area
DiffEqDevTools.best_by_tag
DiffEqDevTools.best_of_families
DiffEqDevTools.with_autodiff_variants
DiffEqDevTools.autoplot
```

## Tableaus

```@docs
DiffEqDevTools.stability_region
DiffEqDevTools.imaginary_stability_interval
DiffEqDevTools.check_tableau
DiffEqDevTools.deduce_Butcher_tableau
DiffEqDevTools.residual_order_condition
```
