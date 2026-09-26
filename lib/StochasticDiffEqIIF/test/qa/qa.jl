using SciMLTesting, StochasticDiffEqIIF, Test
using JET

run_qa(
    StochasticDiffEqIIF;
    reexports_allow = union(public_api_names(StochasticDiffEqCore), (:StochasticDiffEqCore,)),
    explicit_imports = true,
)
