using SciMLTesting, OrdinaryDiffEqFunctionMap, Test

run_qa(
    OrdinaryDiffEqFunctionMap;
    aqua_kwargs = (; piracies = false),  # piracy is needed for default-algorithm dispatch
    explicit_imports = true,
    ei_kwargs = (;
        # OrdinaryDiffEqCore-owned solver-internal hooks that FunctionMap extends but
        # that are not (yet) in OrdinaryDiffEqCore's public-API declaration.
        all_explicit_imports_are_public = (;
            ignore = (:_ode_interpolant!, :isdiscretealg, :isdiscretecache),
        ),
        # SciMLBase-owned default-solve sentinels; non-public in SciMLBase.
        all_qualified_accesses_are_public = (;
            ignore = (:DISCRETE_INPLACE_DEFAULT, :DISCRETE_OUTOFPLACE_DEFAULT),
        ),
    ),
)
