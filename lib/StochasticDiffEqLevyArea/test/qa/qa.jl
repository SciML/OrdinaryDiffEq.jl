using SciMLTesting, StochasticDiffEqLevyArea, Test
using JET

run_qa(
    StochasticDiffEqLevyArea;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(StochasticDiffEqLevyArea), "..", "..", "docs", "src")),
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
)
