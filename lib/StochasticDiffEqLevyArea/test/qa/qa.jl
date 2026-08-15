using SciMLTesting, StochasticDiffEqLevyArea, Test
using JET

run_qa(
    StochasticDiffEqLevyArea;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(StochasticDiffEqLevyArea), "..", "..", "docs", "src")),
    explicit_imports = true,
)
