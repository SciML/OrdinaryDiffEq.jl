using SciMLTesting, OrdinaryDiffEqExplicitTableaus, Test
using JET

run_qa(
    OrdinaryDiffEqExplicitTableaus;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(OrdinaryDiffEqExplicitTableaus), "..", "..", "docs", "src")),
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
)
