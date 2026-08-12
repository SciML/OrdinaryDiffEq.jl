using SciMLTesting, StochasticDiffEqIIF, Test
using JET

run_qa(
    StochasticDiffEqIIF;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(StochasticDiffEqIIF), "..", "..", "docs", "src")),
    reexports_allow = union(public_api_names(StochasticDiffEqCore), (:StochasticDiffEqCore,)),
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
)
