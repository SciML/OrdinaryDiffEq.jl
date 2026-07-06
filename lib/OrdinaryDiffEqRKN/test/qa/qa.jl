using SciMLTesting, OrdinaryDiffEqRKN, Test

run_qa(
    OrdinaryDiffEqRKN;
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # SciMLBase codegen macro, deliberately kept non-public by its owner.
                Symbol("@def"),
                # DiffEqBase perf macro, deliberately kept non-public by its owner.
                Symbol("@tight_loop_macros"),
            ),
        ),
    ),
)
