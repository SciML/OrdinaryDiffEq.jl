using SciMLTesting, StochasticDiffEqLowOrder, Test
using JET

run_qa(
    StochasticDiffEqLowOrder;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(StochasticDiffEqLowOrder), "..", "..", "docs", "src")),
    reexports_allow = union(public_api_names(StochasticDiffEqCore), (:StochasticDiffEqCore,)),
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
    ei_broken = (:no_implicit_imports, :all_explicit_imports_via_owners, :all_explicit_imports_are_public),  # known-broken; see SciML/OrdinaryDiffEq.jl#3776
)
