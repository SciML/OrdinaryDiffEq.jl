using SciMLTesting, OrdinaryDiffEqPRK, Test

run_qa(
    OrdinaryDiffEqPRK;
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore: deliberately non-public codegen macro
                # (kept internal alongside @fold/@OnDemandTableauExtract/@swap!).
                Symbol("@threaded"),
            ),
        ),
    ),
)
