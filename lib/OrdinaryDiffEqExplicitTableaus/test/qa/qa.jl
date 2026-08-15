using SciMLTesting, OrdinaryDiffEqExplicitTableaus, Test
using JET

run_qa(
    OrdinaryDiffEqExplicitTableaus;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(OrdinaryDiffEqExplicitTableaus), "..", "..", "docs", "src")),
    explicit_imports = true,
)
