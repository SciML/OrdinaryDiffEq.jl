using SciMLTesting, DiffEqDevTools, Test

run_qa(
    DiffEqDevTools;
    aqua_kwargs = (;
        ambiguities = false, piracies = false, unbound_args = false, stale_deps = false,
        deps_compat = (; check_extras = false),
    ),
    explicit_imports = true,
    # residual_order_condition is owned by RootedTrees; DiffEqDevTools reexports it
    # after extending it for ODERKTableau (tableau_info.jl).
    reexports_allow = (:residual_order_condition,),
    ei_kwargs = (;
        # SciMLBase-owned solver-interface predicates accessed via SciMLBase (their
        # owner) but not declared `public` in the registered SciMLBase release.
        all_qualified_accesses_are_public = (;
            ignore = (:allowedkeywords, :calculate_ensemble_errors),
        ),
        # SciMLBase-owned `@def`, imported from SciMLBase (its owner) but not declared
        # `public` in the registered SciMLBase release.
        all_explicit_imports_are_public = (;
            ignore = (Symbol("@def"),),
        ),
    ),
)
