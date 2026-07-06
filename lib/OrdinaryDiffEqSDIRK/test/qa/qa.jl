using SciMLTesting, OrdinaryDiffEqSDIRK, Test

# The base solver-author API of OrdinaryDiffEqCore / OrdinaryDiffEqDifferentiation /
# OrdinaryDiffEqNonlinearSolve / DiffEqBase is now declared `public`, so those names
# no longer need ignoring here. What remains are genuine non-public internals with no
# public alternative, grouped by owning package:
#   * OrdinaryDiffEqCore internals that were intentionally NOT promoted to public
#     (the `lorenz`/`lorenz_oop` precompile-workload fixtures and the
#     `trivial_limiter!` default limiter).
#   * SciMLBase private helpers (`_reshape`/`_unwrap_val`/`_vec`).
#   * External packages whose names have no public export (ConstructionBase's
#     `constructorof`, TruncatedStacktraces' `@truncate_stacktrace`).
run_qa(
    OrdinaryDiffEqSDIRK;
    explicit_imports = true,
    ei_kwargs = (
        # `constructorof` is owned by ConstructionBase (not a direct dep) and is
        # reached through SciMLBase's re-export; it is the only non-owner access.
        all_qualified_accesses_via_owners = (; ignore = (:constructorof,)),
        all_qualified_accesses_are_public = (;
            ignore = (
                # non-public in ConstructionBase (via SciMLBase re-export)
                :constructorof,
                # non-public OrdinaryDiffEqCore precompile-workload fixtures
                :lorenz, :lorenz_oop,
            ),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                # non-public SciMLBase internals
                :_reshape, :_unwrap_val, :_vec,
                # non-public TruncatedStacktraces macro
                Symbol("@truncate_stacktrace"),
                # non-public OrdinaryDiffEqCore internal (limiter default)
                :trivial_limiter!,
            ),
        ),
    ),
)
