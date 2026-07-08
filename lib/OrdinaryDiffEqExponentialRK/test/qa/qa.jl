using SciMLTesting, OrdinaryDiffEqExponentialRK, Test

run_qa(
    OrdinaryDiffEqExponentialRK;
    explicit_imports = true,
    ei_kwargs = (
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore owner-internal helpers (deliberately not in the
                # solver-author public API surface declared upstream).
                :_fixup_ad, :fsal_typeof,
                # SciMLBase internal; pending a `public` declaration (SciMLBase#1412).
                :_unwrap_val,
            ),
        ),
    ),
)
