using SciMLTesting, OrdinaryDiffEqImplicitTableaus, Test
using JET

run_qa(
    OrdinaryDiffEqImplicitTableaus;
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
    ei_broken = (:all_qualified_accesses_are_public,),  # known-broken; see SciML/OrdinaryDiffEq.jl#3776
)
