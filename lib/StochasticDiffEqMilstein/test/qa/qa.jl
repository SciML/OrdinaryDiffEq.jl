using SciMLTesting, StochasticDiffEqMilstein, Test
using JET

run_qa(
    StochasticDiffEqMilstein;
    reexports_allow = union(public_api_names(StochasticDiffEqCore), (:StochasticDiffEqCore,)),
    explicit_imports = true,
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            ignore = (
                # StochasticDiffEqLevyArea: Base.ispublic(..., :norv) === false and
                # :norv ∉ names(...); exported equivalent is coefficient_length.
                :norv,
            ),
        ),
    ),
)
