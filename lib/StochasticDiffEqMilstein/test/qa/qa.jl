using SciMLTesting, StochasticDiffEqMilstein, Test
using JET

run_qa(
    StochasticDiffEqMilstein;
    reexports_allow = union(public_api_names(StochasticDiffEqCore), (:StochasticDiffEqCore,)),
    explicit_imports = true,
)
