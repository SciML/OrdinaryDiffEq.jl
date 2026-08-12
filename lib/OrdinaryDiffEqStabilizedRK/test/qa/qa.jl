using SciMLTesting, OrdinaryDiffEqStabilizedRK, Test

run_qa(
    OrdinaryDiffEqStabilizedRK;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(OrdinaryDiffEqStabilizedRK), "..", "..", "docs", "src")),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (
        all_explicit_imports_are_public = (
            ignore = (
                # SciMLBase-owned helpers not yet declared public (pending release).
                :_vec, :value,
            ),
        ),
    ),
)
