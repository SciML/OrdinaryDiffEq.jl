using SciMLTesting, OrdinaryDiffEqImplicitTableaus, Test
using JET

run_qa(
    OrdinaryDiffEqImplicitTableaus;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(OrdinaryDiffEqImplicitTableaus), "..", "..", "docs", "src")),
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
)
