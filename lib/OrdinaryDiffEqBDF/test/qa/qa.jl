using SciMLTesting, OrdinaryDiffEqBDF, Test

run_qa(
    OrdinaryDiffEqBDF;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(OrdinaryDiffEqBDF), "..", "..", "docs", "src")),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore — owner-internal, no public alternative
                :CompositeControllerCache, :lorenz_pref, :lorenz_pref_params,
                # OrdinaryDiffEqNonlinearSolve — owner-internal, no public alternative
                :build_nlsolver, :du_alias_or_new, :markfirststage!, :nlsolve!, :nlsolvefail,
                # Precompile-workload test problems in OrdinaryDiffEqCore, kept
                # non-public on purpose (test fixtures, not solver-author API).
                :lorenz, :lorenz_oop, :lorenz_p, :lorenz_p_params,
            ),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore — owner-internal, no public alternative
                :CompositeControllerCache, :lorenz_pref, :lorenz_pref_params,
                # OrdinaryDiffEqNonlinearSolve — owner-internal, no public alternative
                :build_nlsolver, :du_alias_or_new, :markfirststage!, :nlsolve!, :nlsolvefail,
                # OrdinaryDiffEqCore owner-internal helpers (deliberately excluded
                # from the public solver-author extension surface).
                :_resolved_QT, :trivial_limiter!,
                # SciMLBase internal.
                :_unwrap_val,
                # TruncatedStacktraces internal macro.
                Symbol("@truncate_stacktrace"),
            ),
        ),
    ),
)
