using SciMLTesting, OrdinaryDiffEqFunctionMap, Test

# SciMLBase names re-exported for ordinary ODE usage; everything else stays behind `SciMLBase.`.
const SCIMLBASE_REEXPORTS = (:DiscreteProblem, :solve)

run_qa(
    OrdinaryDiffEqFunctionMap;
    reexports_allow = SCIMLBASE_REEXPORTS,
    aqua_kwargs = (; piracies = false),  # piracy is needed for default-algorithm dispatch
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            # OrdinaryDiffEqCore owner-internal no-op limiter (deliberately not public).
            ignore = (:trivial_limiter!,),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                # SciMLBase-owned default-solve sentinels; non-public in SciMLBase.
                :DISCRETE_INPLACE_DEFAULT,
                :DISCRETE_OUTOFPLACE_DEFAULT,
                # Preserves the statically known NamedTuple type after field removal.
                :structdiff,
            ),
        ),
    ),
)
