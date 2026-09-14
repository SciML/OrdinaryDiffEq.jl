using SciMLTesting, StochasticDiffEqLevyArea, Test
using JET

run_qa(
    StochasticDiffEqLevyArea;
    explicit_imports = true,
)
