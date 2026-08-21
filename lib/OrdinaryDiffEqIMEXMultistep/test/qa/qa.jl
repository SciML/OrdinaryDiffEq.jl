using SciMLTesting, OrdinaryDiffEqIMEXMultistep, SciMLBase, Test

run_qa(
    OrdinaryDiffEqIMEXMultistep;
    # Approve the SciMLBase names this package re-exports. The list itself and the rule
    # behind it are checked repo-wide by test/qa/qa_tests.jl against docs/src/api/reexports.md.
    reexports_allow = intersect(names(SciMLBase), names(OrdinaryDiffEqIMEXMultistep)),
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore: private algorithm-construction helper.
                :_fixup_ad, :_unwrap_val,
                # OrdinaryDiffEqNonlinearSolve owner-internal cross-sublibrary hooks.
                :build_nlsolver, :du_alias_or_new, :markfirststage!, :nlsolve!,
                :nlsolvefail,
            ),
        ),
    ),
)
