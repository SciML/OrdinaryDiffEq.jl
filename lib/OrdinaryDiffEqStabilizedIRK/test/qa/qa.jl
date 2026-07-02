using SciMLTesting, OrdinaryDiffEqStabilizedIRK, Test

run_qa(
    OrdinaryDiffEqStabilizedIRK;
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # SciMLBase internals, not public on registered releases.
                :_unwrap_val, :_reshape, :_vec,
            ),
        ),
    ),
)
