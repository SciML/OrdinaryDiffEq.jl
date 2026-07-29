using SciMLTesting, StochasticDiffEqLevyArea, Test
using JET

run_qa(
    StochasticDiffEqLevyArea;
    # No docs/ tree here; the umbrella manual renders this package's API.
    api_docs_kwargs = (; rendered = false),
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
)
