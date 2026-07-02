using SciMLTesting, ImplicitDiscreteSolve, Test

run_qa(
    ImplicitDiscreteSolve;
    aqua_kwargs = (; piracies = false),
    explicit_imports = true,
    ei_kwargs = (;
        # OrdinaryDiffEqCore solver-internal predicates/hooks that remain
        # non-public: DAE-init entry point, discrete-cache/algorithm predicates,
        # and null-u0 support flag. Not part of the make-public extension surface.
        all_explicit_imports_are_public = (;
            ignore = (
                :isdiscretecache, :isdiscretealg, :_initialize_dae!, :allows_null_u0,
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore controller-resolution internals (non-public).
                :_resolved_QT, :resolve_basic,
                # NonlinearSolveBase solver internals (external, not public).
                :get_fu, :get_u, :not_terminated,
                :update_from_termination_cache!, :update_trace!,
                # SciMLBase initialization-problem predicate (external, not public).
                :has_initializeprob,
            ),
        ),
    ),
)
