using SciMLTesting, OrdinaryDiffEqQPRK, Test

run_qa(
    OrdinaryDiffEqQPRK;
    explicit_imports = true,
    ei_kwargs = (
        all_explicit_imports_are_public = (
            ignore = (
                # OrdinaryDiffEqCore owner-internal codegen macros and default limiter
                # (deliberately not declared public: private codegen/perf surface).
                Symbol("@fold"), Symbol("@OnDemandTableauExtract"), :trivial_limiter!,
            ),
        ),
    ),
)
