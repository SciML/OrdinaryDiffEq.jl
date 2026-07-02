using SciMLTesting, OrdinaryDiffEqPDIRK, Test

run_qa(
    OrdinaryDiffEqPDIRK;
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore: solver-author API is public, but these three
                # are deliberately kept non-public upstream. `@threaded` is a
                # codegen/perf macro; `_fixup_ad` and `differentiation_rk_docstring`
                # are AD-fixup / docstring helpers not part of the extension surface.
                Symbol("@threaded"), :_fixup_ad, :differentiation_rk_docstring,
                # SciMLBase: not declared public on the registered release.
                :_unwrap_val,
            ),
        ),
    ),
)
