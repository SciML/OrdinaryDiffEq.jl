using SciMLTesting, OrdinaryDiffEqRKN, Test

run_qa(
    OrdinaryDiffEqRKN;
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            # Base — owner-internal, no public alternative
            ignore = (:broadcastable,),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                # Base — owner-internal, no public alternative
                :broadcastable,
                # SciMLBase codegen macro, deliberately kept non-public by its owner.
                Symbol("@def"),
                # DiffEqBase perf macro, deliberately kept non-public by its owner.
                Symbol("@tight_loop_macros"),
            ),
        ),
    ),
)
