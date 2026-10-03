# Solver compatibility

- Keep Reactant support in the ordinary solver, controller, and initial-step paths. Use shared traceable control flow; reserve backend checks for tracing boundaries and host-only diagnostics.
- Keep temporary estimates local to traced branches with `let` scopes; avoid exposing logging-macro temporaries as branch outputs.
- Implement dependency-owned methods in the owning package. In particular, specialization policy for SciMLBase function types belongs in SciMLBase.
- `get_fsalfirstlast` initializes storage for composite/default caches; their active FSAL buffers are held by the integrator after algorithm selection. Preserve this distinction when updating derivatives.
- Compiled solves must return the host's `ReturnCode`. Error checks go through `DiffEqBase.staged_check_error`, which evaluates the same predicates as the host's `de_check_error`; never write a second approximation of a check.
- Traced numbers differ from host floats in ways that change failure behavior: `eps(::TracedRNumber)` is `eps` of the type (use `value_eps` / `dt_below_time_eps`), `TracedRNumber` is not `<: Real` so `FastPower.fastpower` falls back to `x^y` (use `controller_fastpower`), and compiled CPU code flushes subnormals to zero. Compare host and compiled retcodes on NaN, overflow, and tiny-step problems when touching these paths.
