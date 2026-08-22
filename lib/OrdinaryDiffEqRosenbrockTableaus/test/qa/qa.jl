using SciMLTesting, OrdinaryDiffEqRosenbrockTableaus, Test
using JET

run_qa(
    OrdinaryDiffEqRosenbrockTableaus;
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
)
