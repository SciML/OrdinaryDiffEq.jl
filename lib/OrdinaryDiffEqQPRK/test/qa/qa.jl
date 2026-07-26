using SciMLTesting, OrdinaryDiffEqQPRK, Test

run_qa(
    OrdinaryDiffEqQPRK;
    # No docs/ tree here; the umbrella manual renders this package's API.
    api_docs_kwargs = (; rendered = false),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (
        all_explicit_imports_are_public = (
            ignore = (
                # OrdinaryDiffEqCore owner-internal codegen macros and default limiter
                # (deliberately not declared public: private codegen/perf surface).
                Symbol("@fold"), Symbol("@OnDemandTableauExtract"), :trivial_limiter!,
            ),
        ),
    ),
)
