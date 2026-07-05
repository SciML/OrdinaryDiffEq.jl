using SciMLTesting, OrdinaryDiffEqNonlinearSolve, Test

run_qa(
    OrdinaryDiffEqNonlinearSolve;
    aqua_kwargs = (; piracies = false),
    explicit_imports = true,
    ei_kwargs = (;
        # Names imported from a re-exporter rather than their defining package.
        # `@SciMLMessage` is owned by SciMLLogging but reached through OrdinaryDiffEqCore
        # (SciMLLogging is not a direct dependency); `WOperator`/`StaticWOperator` are the
        # OrdinaryDiffEqDifferentiation W-operator types, which ExplicitImports attributes
        # to SciMLOperators (their abstract supertype's package).
        all_explicit_imports_via_owners = (;
            ignore = (Symbol("@SciMLMessage"), :WOperator, :StaticWOperator),
        ),
        # Non-public names of upstream packages with no public alternative. Each is owned
        # by the package it is imported from; remove an entry once that package marks the
        # name `public`/exports it. Tracked in SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                # `@SciMLMessage` reached through OrdinaryDiffEqCore (owner SciMLLogging).
                Symbol("@SciMLMessage"),
                # SciMLBase internals (owner of these names but not public).
                :_vec, :_reshape, :postamble!,
                # SciMLOperators abstract type (not public).
                :AbstractSciMLOperator,
                # OrdinaryDiffEqDifferentiation W-operator types (attributed to SciMLOperators).
                :WOperator, :StaticWOperator,
                # ForwardDiff / StaticArraysCore internals.
                :Dual, :StaticArray,
            ),
        ),
        # Qualified `Module.name` accesses to non-public names with no public alternative.
        all_qualified_accesses_are_public = (;
            ignore = (
                # SciMLBase internals.
                :value, :anyeltypedual, :forwarddiff_chunksize,
                :has_Wfact, :has_Wfact_t, :has_jac_u, :has_jac_du, :postamble!,
                # ForwardDiff internals.
                :Tag, :pickchunksize,
            ),
        ),
    ),
)
