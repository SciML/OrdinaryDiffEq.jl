using SciMLTesting, OrdinaryDiffEqImplicitTableaus, Test
using JET

run_qa(
    OrdinaryDiffEqImplicitTableaus;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(OrdinaryDiffEqImplicitTableaus), "..", "..", "docs", "src")),
    explicit_imports = true,
)
