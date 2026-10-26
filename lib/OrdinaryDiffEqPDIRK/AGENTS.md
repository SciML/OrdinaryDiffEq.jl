# PDIRK concurrency

- Nonlinear workers own their statistics buffers and Jacobian evaluation state. Do not write shared integrator statistics or `force_stepfail` from workers, or add a global lock.
- After each threaded loop joins, the parent merges buffers in solver order and assigns the aggregate failure flag before any early return.
- Serial stepping uses one nonlinear solver for all four stages; threaded stepping uses two solvers and merges their private buffers after each join. Retained out-of-place stage values must not alias storage reused by a subsequent nonlinear solve.
- Concurrency changes must pass the package's `threaded_stats_tests.jl` through the existing test runner with at least four Julia threads. Cover supported nonlinear solvers, every statistics field, and failure paths; demonstrate the regression on the unfixed code.
