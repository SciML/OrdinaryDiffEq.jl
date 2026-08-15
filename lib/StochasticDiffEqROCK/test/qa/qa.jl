using SciMLTesting, StochasticDiffEqROCK, Test
using JET

run_qa(
    StochasticDiffEqROCK;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(StochasticDiffEqROCK), "..", "..", "docs", "src")),
    reexports_allow = union(public_api_names(StochasticDiffEqCore), (:StochasticDiffEqCore,)),
    explicit_imports = true,
    ei_broken = (:no_implicit_imports, :no_stale_explicit_imports, :all_explicit_imports_via_owners, :all_explicit_imports_are_public),  # known-broken; see SciML/OrdinaryDiffEq.jl#3776
)
