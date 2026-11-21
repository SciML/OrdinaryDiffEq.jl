using SciMLTesting, OrdinaryDiffEqSymplecticRK, Test

run_qa(
    OrdinaryDiffEqSymplecticRK;
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    # the tstop a step lands on (as in DelayDiffEq's QA ignores)
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            ignore = (:_get_next_step_tstop, :_get_tstop_target),
        ),
    ),
)
