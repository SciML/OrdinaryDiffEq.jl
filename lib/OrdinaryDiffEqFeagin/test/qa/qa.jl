using SciMLTesting, OrdinaryDiffEqFeagin, Test

# SciMLBase names re-exported for ordinary ODE usage; everything else stays behind `SciMLBase.`.
const SCIMLBASE_REEXPORTS = (
    :ODEProblem, :ODEFunction, :solve, :init, :solve!, :step!, :remake, :reinit!,
    :ReturnCode, :ContinuousCallback, :DiscreteCallback, :VectorContinuousCallback,
    :CallbackSet, :terminate!,
)

run_qa(
    OrdinaryDiffEqFeagin;
    reexports_allow = SCIMLBASE_REEXPORTS,
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
