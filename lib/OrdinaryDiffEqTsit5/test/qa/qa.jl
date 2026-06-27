using SciMLTesting, OrdinaryDiffEqTsit5, Test

run_qa(
    OrdinaryDiffEqTsit5;
    explicit_imports = true,
    ei_broken = (:no_implicit_imports, :all_explicit_imports_via_owners, :all_qualified_accesses_are_public, :all_explicit_imports_are_public),  # known-broken; see SciML/OrdinaryDiffEq.jl#3776
)
