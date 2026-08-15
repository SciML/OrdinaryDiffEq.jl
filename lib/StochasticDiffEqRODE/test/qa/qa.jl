using SciMLTesting, StochasticDiffEqRODE, Test
using JET

run_qa(
    StochasticDiffEqRODE;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(StochasticDiffEqRODE), "..", "..", "docs", "src")),
    reexports_allow = union(public_api_names(StochasticDiffEqCore), (:StochasticDiffEqCore,)),
    explicit_imports = true,
    ei_broken = (:no_implicit_imports, :all_explicit_imports_via_owners, :all_explicit_imports_are_public),  # known-broken; see SciML/OrdinaryDiffEq.jl#3776
)
