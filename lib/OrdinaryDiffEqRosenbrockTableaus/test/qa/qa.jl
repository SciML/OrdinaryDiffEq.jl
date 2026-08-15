using SciMLTesting, OrdinaryDiffEqRosenbrockTableaus, Test
using JET

run_qa(
    OrdinaryDiffEqRosenbrockTableaus;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(OrdinaryDiffEqRosenbrockTableaus), "..", "..", "docs", "src")),
    explicit_imports = true,
)
