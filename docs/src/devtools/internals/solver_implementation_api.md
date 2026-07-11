# Solver Implementation Extension API

These extension points support the solver sublibraries in the OrdinaryDiffEq
monorepo and external solver implementations.

## Nonlinear solver hooks

The nonlinear-solver implementation package provides the construction, iteration,
status, and fixed-point hooks used by implicit solver sublibraries.

```@docs
OrdinaryDiffEqCore.AbstractNLSolverCache
OrdinaryDiffEqNonlinearSolve.build_nlsolver
OrdinaryDiffEqNonlinearSolve.nlsolve!
OrdinaryDiffEqNonlinearSolve.nlsolvefail
OrdinaryDiffEqNonlinearSolve.compute_step!
OrdinaryDiffEqNonlinearSolve.initial_η
OrdinaryDiffEqNonlinearSolve.markfirststage!
OrdinaryDiffEqNonlinearSolve.du_alias_or_new
OrdinaryDiffEqNonlinearSolve.anderson
OrdinaryDiffEqNonlinearSolve.anderson!
```

## Solver caches

Solver sublibraries define their cache types and construct them through these
extension points. Automatic and composite algorithms use the concrete cache types
listed here to delegate between their constituent algorithms.

```@docs
OrdinaryDiffEqCore.OrdinaryDiffEqConstantCache
OrdinaryDiffEqCore.OrdinaryDiffEqMutableCache
OrdinaryDiffEqCore.DefaultCache
OrdinaryDiffEqCore.AutoSwitchCache
OrdinaryDiffEqCore.CompositeControllerCache
OrdinaryDiffEqCore.PredictiveControllerCache
OrdinaryDiffEqCore.alg_cache
OrdinaryDiffEqCore.get_fsalfirstlast
OrdinaryDiffEqCore.handle_callback_modifiers!
OrdinaryDiffEqCore.@cache
```

## Jacobian / W-matrix / differentiation configuration

Provided by `OrdinaryDiffEqDifferentiation`: the Jacobian, W-matrix,
time-derivative, and linear-solve hooks that implicit solver sublibraries build
their caches and `perform_step!` implementations on.

```@docs
OrdinaryDiffEqDifferentiation.build_J_W
OrdinaryDiffEqDifferentiation.build_uf
OrdinaryDiffEqDifferentiation.build_jac_config
OrdinaryDiffEqDifferentiation.build_grad_config
OrdinaryDiffEqDifferentiation.calc_J
OrdinaryDiffEqDifferentiation.calc_J!
OrdinaryDiffEqDifferentiation.calc_tderivative
OrdinaryDiffEqDifferentiation.calc_tderivative!
OrdinaryDiffEqDifferentiation.calc_rosenbrock_differentiation
OrdinaryDiffEqDifferentiation.calc_rosenbrock_differentiation!
OrdinaryDiffEqDifferentiation.jacobian!
OrdinaryDiffEqDifferentiation.jacobian2W!
OrdinaryDiffEqDifferentiation.update_W!
OrdinaryDiffEqDifferentiation.resize_jac_config!
OrdinaryDiffEqDifferentiation.resize_grad_config!
OrdinaryDiffEqDifferentiation.dolinsolve
OrdinaryDiffEqDifferentiation.wrapprecs
OrdinaryDiffEqDifferentiation.is_always_new
OrdinaryDiffEqDifferentiation.islinearfunction
OrdinaryDiffEqDifferentiation.issuccess_W
```
