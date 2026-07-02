using SciMLTesting, StochasticDiffEqIIF, Test
using JET

run_qa(
    StochasticDiffEqIIF;
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
)
