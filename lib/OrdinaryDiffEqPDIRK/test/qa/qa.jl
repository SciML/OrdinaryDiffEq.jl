using SciMLTesting, OrdinaryDiffEqPDIRK, Test

# `public` on a name another package owns counts as a public reexport to
# SciMLTesting, so the threading options need approving here too.
const THREADING_PUBLIC = (:Sequential, :BaseThreads, :PolyesterThreads)

# SciMLBase names re-exported for ordinary ODE usage; everything else stays behind `SciMLBase.`.
const SCIMLBASE_REEXPORTS = (
    :ODEProblem, :ODEFunction, :solve, :init, :solve!, :step!, :remake, :reinit!,
    :ReturnCode, :ContinuousCallback, :DiscreteCallback, :VectorContinuousCallback,
    :CallbackSet, :terminate!,
)

run_qa(
    OrdinaryDiffEqPDIRK;
    reexports_allow = (SCIMLBASE_REEXPORTS..., THREADING_PUBLIC...),
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
