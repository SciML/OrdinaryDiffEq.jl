using SciMLTesting, OrdinaryDiffEqExplicitRK, Test

run_qa(
    OrdinaryDiffEqExplicitRK;
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # TruncatedStacktraces-owned macro (external, not `public`).
                Symbol("@truncate_stacktrace"),
                # OrdinaryDiffEqCore owner-internal no-op limiter (deliberately not public).
                :trivial_limiter!,
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                # Base.Cartesian codegen macros (external, not `public`).
                Symbol("@nexprs"), Symbol("@nif"),
            ),
        ),
    ),
    api_docs_kwargs = (;
        # Reexported upstream SciMLOperators names are documented at their owner.
        ignore = (:StaticWOperator, :has_concretization),
    ),
)
