using SciMLTesting, OrdinaryDiffEqDefault, Test

run_qa(
    OrdinaryDiffEqDefault;
    aqua_kwargs = (; piracies = false),
    explicit_imports = true,
    ei_broken = (:no_implicit_imports, :no_stale_explicit_imports, :all_qualified_accesses_are_public, :all_explicit_imports_are_public),  # known-broken; see SciML/OrdinaryDiffEq.jl#3776
)
