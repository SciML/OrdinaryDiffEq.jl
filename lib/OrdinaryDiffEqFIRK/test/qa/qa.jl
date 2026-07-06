using SciMLTesting, OrdinaryDiffEqFIRK, Test

# The solver-author API of OrdinaryDiffEqCore / OrdinaryDiffEqDifferentiation /
# OrdinaryDiffEqNonlinearSolve / DiffEqBase is now declared `public` on this branch,
# so those names no longer need ignoring. The residual below is the genuine set of
# non-public internals with no public alternative, grouped by owning package. Any
# name later made public will surface as an Unexpected Pass and can be dropped.
# See SciML/OrdinaryDiffEq.jl#3776.
run_qa(
    OrdinaryDiffEqFIRK;
    explicit_imports = true,
    ei_kwargs = (
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore — private codegen macro / default no-op limiter,
                # kept owner-internal (no public alternative).
                Symbol("@threaded"), :trivial_limiter!,
                # SciMLBase internal helpers (pending public decl, SciMLBase#1412).
                :_reshape, :_unwrap_val, :_vec, :value,
                # Genuine external deps, non-public in their owner.
                :fastpower,               # FastPower
                :AbstractSciMLOperator,   # SciMLOperators
            ),
        ),
    ),
)
