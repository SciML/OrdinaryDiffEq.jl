using SciMLTesting, StochasticDiffEqLevyArea, Test
using JET

run_qa(
    StochasticDiffEqLevyArea;
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
)
