using SciMLTesting, OrdinaryDiffEqRosenbrockTableaus, Test
using JET

run_qa(
    OrdinaryDiffEqRosenbrockTableaus;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(OrdinaryDiffEqRosenbrockTableaus), "..", "..", "docs", "src")),
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
)
