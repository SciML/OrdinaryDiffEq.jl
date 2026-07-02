using SciMLTesting, OrdinaryDiffEqBDF, Test

run_qa(
    OrdinaryDiffEqBDF;
    explicit_imports = true,
    ei_kwargs = (
        all_qualified_accesses_are_public = (;
            ignore = (
                # Precompile-workload test problems in OrdinaryDiffEqCore, kept
                # non-public on purpose (test fixtures, not solver-author API).
                :lorenz, :lorenz_oop,
            ),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
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
