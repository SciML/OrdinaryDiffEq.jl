using SciMLTesting, OrdinaryDiffEqAdamsBashforthMoulton, Test

run_qa(
    OrdinaryDiffEqAdamsBashforthMoulton;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(OrdinaryDiffEqAdamsBashforthMoulton), "..", "..", "docs", "src")),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internal limiter hook (owner-internal;
                # deliberately not declared public, like the @fold/@threaded codegen macros).
                :trivial_limiter!,
            ),
        ),
    ),
)
