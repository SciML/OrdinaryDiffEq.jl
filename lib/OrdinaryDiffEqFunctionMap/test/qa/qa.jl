using SciMLTesting, OrdinaryDiffEqFunctionMap, Test

run_qa(
    OrdinaryDiffEqFunctionMap;
    aqua_kwargs = (; piracies = false),  # piracy is needed for default-algorithm dispatch
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            # OrdinaryDiffEqCore owner-internal no-op limiter (deliberately not public).
            ignore = (:trivial_limiter!,),
        ),
        # SciMLBase-owned default-solve sentinels; non-public in SciMLBase.
        all_qualified_accesses_are_public = (;
            ignore = (:DISCRETE_INPLACE_DEFAULT, :DISCRETE_OUTOFPLACE_DEFAULT),
        ),
    ),
)
