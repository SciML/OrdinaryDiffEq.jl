using SciMLTesting, OrdinaryDiffEqAMF, SciMLBase, Test

run_qa(
    OrdinaryDiffEqAMF;
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            # LinearSolve does not declare these extension points and cache internals public.
            ignore = (
                :LinearCache, :SciMLLinearSolveAlgorithm, :init_cacheval, :needs_concrete_A,
            ),
        ),
    ),
)
