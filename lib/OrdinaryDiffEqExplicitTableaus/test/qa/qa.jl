using SciMLTesting, OrdinaryDiffEqExplicitTableaus, Test
using JET

run_qa(
    OrdinaryDiffEqExplicitTableaus;
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
)
