using SciMLTesting, OrdinaryDiffEqQPRK, SciMLBase, Test

run_qa(
    OrdinaryDiffEqQPRK;
    # Approve the SciMLBase names this package re-exports. The list itself and the rule
    # behind it are checked repo-wide by test/qa/qa_tests.jl against docs/src/api/reexports.md.
    reexports_allow = intersect(names(SciMLBase), names(OrdinaryDiffEqQPRK)),
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
