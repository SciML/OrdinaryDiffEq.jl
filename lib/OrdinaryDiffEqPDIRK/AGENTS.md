# PDIRK concurrency

- Nonlinear workers own their statistics buffers and Jacobian evaluation state. Do not write shared integrator statistics or `force_stepfail` from workers, or add a global lock.
- After each threaded loop joins, the parent merges buffers in solver order and assigns the aggregate failure flag before any early return.
- Preserve the same two-solver stage ordering in serial and threaded modes. Retained out-of-place stage values must not alias storage reused by a subsequent nonlinear solve.
- Concurrency changes must pass the package's `threaded_stats_tests.jl` through the existing test runner with at least four Julia threads. Cover supported nonlinear solvers, every statistics field, and failure paths; demonstrate the regression on the unfixed code.
