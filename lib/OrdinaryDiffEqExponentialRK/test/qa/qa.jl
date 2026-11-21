using SciMLTesting, OrdinaryDiffEqExponentialRK, Test

run_qa(
    OrdinaryDiffEqExponentialRK;
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (
        all_qualified_accesses_are_public = (;
            ignore = (
                # Base — owner-internal, no public alternative
                :structdiff,
                # OrdinaryDiffEqCore — internal dispatch hook extended for ETD2 caches
                :reset_fsal!,
            ),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                # Base — owner-internal, no public alternative
                :structdiff,
                # OrdinaryDiffEqCore owner-internal helpers (deliberately not in the
                # solver-author public API surface declared upstream).
                :_fixup_ad, :fsal_typeof, :trivial_limiter!,
                # SciMLBase internal; pending a `public` declaration (SciMLBase#1412).
                :_unwrap_val,
            ),
        ),
    ),
)
