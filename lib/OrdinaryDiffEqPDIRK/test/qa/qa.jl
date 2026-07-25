using SciMLTesting, OrdinaryDiffEqPDIRK, Test

run_qa(
    OrdinaryDiffEqPDIRK;
    # No docs/ tree here; the umbrella manual renders this package's API.
    api_docs_kwargs = (; rendered = false),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore: `@threaded` is a private codegen/perf macro,
                # deliberately kept non-public upstream (owner-internal).
                Symbol("@threaded"),
                # OrdinaryDiffEqNonlinearSolve owner-internal cross-sublibrary hooks;
                # no public wrapper exists.
                :build_nlsolver, :markfirststage!, :nlsolve!, :nlsolvefail,
                # SciMLBase: not declared public on the registered release.
                :_unwrap_val,
            ),
        ),
    ),
)
