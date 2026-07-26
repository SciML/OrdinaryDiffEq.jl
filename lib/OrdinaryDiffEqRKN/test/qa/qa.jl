using SciMLTesting, OrdinaryDiffEqRKN, Test

run_qa(
    OrdinaryDiffEqRKN;
    # No docs/ tree here; the umbrella manual renders this package's API.
    api_docs_kwargs = (; rendered = false),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # SciMLBase codegen macro, deliberately kept non-public by its owner.
                Symbol("@def"),
                # DiffEqBase perf macro, deliberately kept non-public by its owner.
                Symbol("@tight_loop_macros"),
            ),
        ),
    ),
)
