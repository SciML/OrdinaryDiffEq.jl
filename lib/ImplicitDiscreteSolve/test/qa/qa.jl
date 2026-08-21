using SciMLTesting, ImplicitDiscreteSolve, Test

run_qa(
    ImplicitDiscreteSolve;
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    aqua_kwargs = (; piracies = false),
    explicit_imports = true,
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore controller-resolution internal (owner-internal,
                # deliberately kept non-public in the make-public extension surface).
                :_resolved_QT,
                # SciMLBase initialization-problem predicate (external, not public).
                :has_initializeprob,
            ),
        ),
    ),
)
