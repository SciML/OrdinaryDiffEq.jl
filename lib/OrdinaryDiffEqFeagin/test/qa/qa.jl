using SciMLTesting, OrdinaryDiffEqFeagin, SciMLBase, Test

run_qa(
    OrdinaryDiffEqFeagin;
    # Approve the SciMLBase names this package re-exports. The list itself and the rule
    # behind it are checked repo-wide by test/qa/qa_tests.jl against docs/src/api/reexports.md.
    reexports_allow = intersect(names(SciMLBase), names(OrdinaryDiffEqFeagin)),
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore-owned internals, deliberately not `public`.
                :CompiledFloats, :trivial_limiter!,
                # DiffEqBase-owned internal macro, deliberately not `public`.
                Symbol("@tight_loop_macros"),
            ),
        ),
    ),
)
