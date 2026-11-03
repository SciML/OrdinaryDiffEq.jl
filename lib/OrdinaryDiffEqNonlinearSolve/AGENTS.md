# Nonlinear solver validation

- Stalled line-search safeguards use the undamped correction in the same tolerance
  units and algorithm-specific error scaling as `compute_step!`. Cover both
  substantial corrections with unchanged residuals and legitimate tiny corrections.
- Relaxation regressions cover cached and fresh Jacobians, assert stage residual
  safety, and include independent large state components.
- Line-search trial RHS calls do not increment `nf`. Justify work budgets with
  measured solver statistics. Agreement between Float64 solves is an agreement
  metric; accuracy claims use an independent, converged reference solution.
- Run the package's Core and QA groups through its existing test harness.
