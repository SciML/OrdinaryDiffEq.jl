using SciMLTesting, OrdinaryDiffEqSymplecticRK, Test

run_qa(
    OrdinaryDiffEqSymplecticRK;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(OrdinaryDiffEqSymplecticRK), "..", "..", "docs", "src")),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
)
