using SciMLTesting, OrdinaryDiffEqStabilizedIRK, Test

run_qa(
    OrdinaryDiffEqStabilizedIRK;
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore owner-internal helpers deliberately kept
                # non-public (perf/dispatch/misc utilities, not part of the
                # documented solver-author extension surface).
                :fac_default_gamma, :_fixup_ad, :isnewton,
                # SciMLBase internals, not public on registered releases.
                :_unwrap_val, :_reshape, :_vec,
            ),
        ),
    ),
)
