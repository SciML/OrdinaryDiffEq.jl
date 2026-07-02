using SciMLTesting, ImplicitDiscreteSolve, Test

run_qa(
    ImplicitDiscreteSolve;
    aqua_kwargs = (; piracies = false),
    explicit_imports = true,
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore controller-resolution internal (owner-internal,
                # deliberately kept non-public in the make-public extension surface).
                :_resolved_QT,
                # NonlinearSolveBase solver internals (external, not public).
                :get_fu, :get_u, :not_terminated,
                :update_from_termination_cache!, :update_trace!,
                # SciMLBase initialization-problem predicate (external, not public).
                :has_initializeprob,
            ),
        ),
    ),
)
