using SciMLTesting, OrdinaryDiffEqQPRK, Test

# SciMLBase names re-exported for ordinary ODE usage; everything else stays behind `SciMLBase.`.
const SCIMLBASE_REEXPORTS = (
    :ODEProblem, :ODEFunction, :solve, :init, :solve!, :step!, :remake, :reinit!,
    :ReturnCode, :ContinuousCallback, :DiscreteCallback, :VectorContinuousCallback,
    :CallbackSet, :terminate!,
)

run_qa(
    OrdinaryDiffEqQPRK;
    reexports_allow = SCIMLBASE_REEXPORTS,
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
