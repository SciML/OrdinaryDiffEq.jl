using SciMLTesting, OrdinaryDiffEqPDIRK, Test

run_qa(
    OrdinaryDiffEqPDIRK;
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore: `@threaded` is a private codegen/perf macro,
                # deliberately kept non-public upstream (owner-internal).
                Symbol("@threaded"),
                # SciMLBase: not declared public on the registered release.
                :_unwrap_val,
            ),
        ),
    ),
)
