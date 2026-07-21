using SciMLTesting, OrdinaryDiffEqIMEXMultistep, Test

run_qa(
    OrdinaryDiffEqIMEXMultistep;
    explicit_imports = true,
    ei_kwargs = (;
        # `@reexport using SciMLBase` brings the `SciMLBase` module name into scope
        # implicitly; this is the intended SciML re-export pattern, not a stray import.
        no_implicit_imports = (; ignore = (:SciMLBase,)),
        # Solver-internal helpers imported from their owning packages but deliberately
        # kept non-public there:
        #   OrdinaryDiffEqCore: `_fixup_ad` (autodiff-fixup private helper)
        #   SciMLBase:          `_unwrap_val`
        # OrdinaryDiffEqNonlinearSolve owner-internal cross-sublibrary hooks;
        # no public wrapper exists.
        all_explicit_imports_are_public = (;
            ignore = (
                :_fixup_ad, :_unwrap_val,
                :build_nlsolver, :du_alias_or_new, :markfirststage!, :nlsolve!,
                :nlsolvefail,
            ),
        ),
    ),
)
