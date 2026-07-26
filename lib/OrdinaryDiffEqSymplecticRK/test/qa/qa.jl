using SciMLTesting, OrdinaryDiffEqSymplecticRK, Test

run_qa(
    OrdinaryDiffEqSymplecticRK;
    # No docs/ tree here; the umbrella manual renders this package's API.
    api_docs_kwargs = (; rendered = false),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
)
