using SciMLTesting, OrdinaryDiffEqIMEXMultistep, Test

# SciMLBase names re-exported for ordinary ODE usage; everything else stays behind `SciMLBase.`.
const SCIMLBASE_REEXPORTS = (:SplitODEProblem, :solve)

run_qa(
    OrdinaryDiffEqIMEXMultistep;
    reexports_allow = SCIMLBASE_REEXPORTS,
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
