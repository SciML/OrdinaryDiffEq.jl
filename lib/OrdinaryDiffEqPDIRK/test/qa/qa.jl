using SciMLTesting, OrdinaryDiffEqPDIRK, SciMLBase, Test

# `public` on a name another package owns counts as a public reexport to
# SciMLTesting, so the threading options need approving here too.
const THREADING_PUBLIC = (:Sequential, :BaseThreads, :PolyesterThreads)

run_qa(
    OrdinaryDiffEqPDIRK;
    # Approve the SciMLBase names this package re-exports. The list itself and the rule
    # behind it are checked repo-wide by test/qa/qa_tests.jl against docs/src/api/reexports.md.
    reexports_allow = vcat(intersect(names(SciMLBase), names(OrdinaryDiffEqPDIRK)), collect(THREADING_PUBLIC)),
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore: `@threaded` is a private codegen/perf macro,
                # deliberately kept non-public upstream (owner-internal).
                Symbol("@threaded"),
                # OrdinaryDiffEqNonlinearSolve owner-internal cross-sublibrary hooks;
                # no public wrapper exists.
                :build_nlsolver, :markfirststage!, :nlsolve!, :nlsolvefail,
                # SciMLBase: not declared public on the registered release.
                :_unwrap_val,
            ),
        ),
    ),
)
