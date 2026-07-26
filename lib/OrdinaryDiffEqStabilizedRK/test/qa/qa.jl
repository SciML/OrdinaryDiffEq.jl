using SciMLTesting, OrdinaryDiffEqStabilizedRK, Test

run_qa(
    OrdinaryDiffEqStabilizedRK;
    # No docs/ tree here; the umbrella manual renders this package's API.
    api_docs_kwargs = (; rendered = false),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (
        all_explicit_imports_are_public = (
            ignore = (
                # SciMLBase-owned helpers not yet declared public (pending release).
                :_vec, :value,
                # OrdinaryDiffEqCore owner-internal no-op limiter (deliberately not public).
                :trivial_limiter!,
            ),
        ),
    ),
)
