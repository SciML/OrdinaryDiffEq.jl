using SciMLTesting, StochasticDiffEqIIF, Test
using JET

run_qa(
    StochasticDiffEqIIF;
    reexports_allow = union(public_api_names(StochasticDiffEqCore), (:StochasticDiffEqCore,)),
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
)
