# Nonlinear convergence

- Use independently derived implicit-stage equations for solver-level globalization regressions. Check the stage residual and solution when a solve succeeds, and verify explicit failure when the stage cannot be accepted.
- A small damped Newton correction alone does not establish convergence. When damping makes a large correction appear converged, independently evaluate the accepted stage residual. Apply this requirement to every convergence exit, including convergence-rate estimates.
- Preserve valid sub-epsilon line-search step lengths. Judge accepted steps by representable iterate changes and residual progress rather than an absolute epsilon cutoff on the step length.
