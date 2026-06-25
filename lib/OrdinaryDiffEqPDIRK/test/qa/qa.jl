using SciMLTesting, OrdinaryDiffEqPDIRK, Test

run_qa(
    OrdinaryDiffEqPDIRK;
    explicit_imports = true,
    ei_broken = (:no_implicit_imports, :no_stale_explicit_imports, :all_explicit_imports_via_owners, :all_explicit_imports_are_public),  # known-broken; see SciML/OrdinaryDiffEq.jl#3776
)
