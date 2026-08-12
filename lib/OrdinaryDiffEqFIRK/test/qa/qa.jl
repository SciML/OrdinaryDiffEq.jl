using SciMLTesting, OrdinaryDiffEqFIRK, Test

# `public` on a name another package owns counts as a public reexport to
# SciMLTesting, so the threading options need approving here too.
const THREADING_PUBLIC = (:Sequential, :BaseThreads, :PolyesterThreads)

# The solver-author API of OrdinaryDiffEqCore / OrdinaryDiffEqDifferentiation /
# OrdinaryDiffEqNonlinearSolve / DiffEqBase is now declared `public` on this branch,
# so those names no longer need ignoring. The residual below is the genuine set of
# non-public internals with no public alternative, grouped by owning package. Any
# name later made public will surface as an Unexpected Pass and can be dropped.
# See SciML/OrdinaryDiffEq.jl#3776.
run_qa(
    OrdinaryDiffEqFIRK;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(OrdinaryDiffEqFIRK), "..", "..", "docs", "src")),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,), THREADING_PUBLIC),
    explicit_imports = true,
    ei_kwargs = (
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore — owner-internal, no public alternative
                :PredictiveControllerCache,
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
        all_qualified_accesses_are_public = (;
            ignore = (
                # LinearSolve internal, and the only way to ask whether a linear solver
                # can consume an operator `A`. `build_J_W` in
                # OrdinaryDiffEqDifferentiation reaches for the same name and ignores it
                # the same way; FIRK has to make the matching decision for its own
                # shifted stage matrices. Make-public candidate upstream.
                :needs_concrete_A,
            ),
        ),
    ),
)
