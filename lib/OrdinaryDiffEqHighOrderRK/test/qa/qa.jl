using SciMLTesting, OrdinaryDiffEqHighOrderRK, Test

run_qa(
    OrdinaryDiffEqHighOrderRK;
    explicit_imports = true,
    ei_kwargs = (
        # `@reexport using SciMLBase` is the intended public-API reexport; it
        # necessarily brings the `SciMLBase`/`Reexport`/`@reexport` bindings in
        # implicitly, which there is no way to make explicit.
        no_implicit_imports = (ignore = (:SciMLBase, :Reexport, Symbol("@reexport")),),
        # `@tight_loop_macros` is used inside the `@muladd`-expanded perform_step/
        # addsteps bodies, which ExplicitImports cannot see, so it reads as a stale
        # positive; it is a DiffEqBase-internal codegen macro (not public upstream).
        no_stale_explicit_imports = (ignore = (Symbol("@tight_loop_macros"),),),
        all_explicit_imports_are_public = (
            ignore = (
                # DiffEqBase-internal codegen macro (not declared public upstream).
                Symbol("@tight_loop_macros"),
                # OrdinaryDiffEqCore-internal no-op limiter default, deliberately
                # kept non-public.
                :trivial_limiter!,
            ),
        ),
    ),
)
