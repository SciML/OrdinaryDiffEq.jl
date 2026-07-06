using SciMLTesting, StochasticDiffEqRODE, Test
using JET

run_qa(
    StochasticDiffEqRODE;
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
    ei_broken = (:no_implicit_imports, :all_explicit_imports_via_owners, :all_qualified_accesses_via_owners, :all_explicit_imports_are_public),  # known-broken; see SciML/OrdinaryDiffEq.jl#3776
)
