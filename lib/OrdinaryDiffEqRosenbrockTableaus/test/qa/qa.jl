using SciMLTesting, OrdinaryDiffEqRosenbrockTableaus, Test
using JET

run_qa(
    OrdinaryDiffEqRosenbrockTableaus;
    # No docs/ tree here; the umbrella manual renders this package's API.
    api_docs_kwargs = (; rendered = false),
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
)
