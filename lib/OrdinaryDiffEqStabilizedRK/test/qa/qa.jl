using SciMLTesting, OrdinaryDiffEqStabilizedRK, Test

run_qa(
    OrdinaryDiffEqStabilizedRK;
    explicit_imports = true,
    ei_kwargs = (
        all_explicit_imports_are_public = (
            ignore = (
                # SciMLBase-owned helpers not yet declared public (pending release).
                :_vec, :value,
            ),
        ),
    ),
)
