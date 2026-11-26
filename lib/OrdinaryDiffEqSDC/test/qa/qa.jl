using SciMLTesting, OrdinaryDiffEqSDC, Test

run_qa(
    OrdinaryDiffEqSDC;
    # No docs/ tree here; the umbrella manual renders this package's API.
    api_docs_kwargs = (; rendered = false),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (
        all_explicit_imports_are_public = (;
            ignore = (Symbol("@threaded"),),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                # non-public OrdinaryDiffEqCore precompile-workload fixtures
                :lorenz, :lorenz_oop, :lorenz_p, :lorenz_p_params,
            ),
        ),
    ),
)
