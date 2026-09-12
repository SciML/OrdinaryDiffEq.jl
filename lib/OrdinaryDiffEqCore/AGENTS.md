# Solver compatibility

- Keep Reactant support in the ordinary solver, controller, and initial-step paths. Use shared traceable control flow; reserve backend checks for tracing boundaries and host-only diagnostics.
- Keep temporary estimates local to traced branches with `let` scopes; avoid exposing logging-macro temporaries as branch outputs.
- Implement dependency-owned methods in the owning package. In particular, specialization policy for SciMLBase function types belongs in SciMLBase.
- `get_fsalfirstlast` initializes storage for composite/default caches; their active FSAL buffers are held by the integrator after algorithm selection. Preserve this distinction when updating derivatives.
