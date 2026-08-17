using SciMLTesting, OrdinaryDiffEqImplicitTableaus, Test
using JET

run_qa(
    OrdinaryDiffEqImplicitTableaus;
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
)
