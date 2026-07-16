using SciMLTesting, OrdinaryDiffEqFeagin, Test

run_qa(
    OrdinaryDiffEqFeagin;
    explicit_imports = true,
    ei_kwargs = (;
        # `@reexport using SciMLBase` deliberately re-exports SciMLBase's API
        # as part of this solver package's public surface.
        no_implicit_imports = (; ignore = (:SciMLBase,)),
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
