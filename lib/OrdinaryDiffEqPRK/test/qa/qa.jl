using SciMLTesting, OrdinaryDiffEqPRK, Test

run_qa(
    OrdinaryDiffEqPRK;
    # No docs/ tree here; the umbrella manual renders this package's API.
    api_docs_kwargs = (; rendered = false),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore: deliberately non-public codegen macro
                # (kept internal alongside @fold/@OnDemandTableauExtract/@swap!).
                Symbol("@threaded"),
            ),
        ),
    ),
)
