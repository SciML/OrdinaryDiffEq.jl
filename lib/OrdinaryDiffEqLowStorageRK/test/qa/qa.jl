using SciMLTesting, OrdinaryDiffEqLowStorageRK, Test

run_qa(
    OrdinaryDiffEqLowStorageRK;
    # No docs/ tree here; the umbrella manual renders this package's API.
    api_docs_kwargs = (; rendered = false),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (
        # Every remaining name is a genuine non-public internal of its owner. The
        # solver-author API of OrdinaryDiffEqCore / DiffEqBase is now declared
        # `public`, so those ignores were dropped.
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore — owner-internal, no public alternative
                :lorenz_pref, :lorenz_pref_params,
                # OrdinaryDiffEqCore precompile-workload probes (owner-internal)
                :lorenz, :lorenz_oop, :lorenz_p, :lorenz_p_params,
                # Base.Broadcast internals used in ArrayFuse copyto!/materialize!
                :Broadcasted, :materialize!,
            ),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore — owner-internal, no public alternative
                :lorenz_pref, :lorenz_pref_params,
                # OrdinaryDiffEqCore default no-op limiter (owner-internal)
                :trivial_limiter!,
            ),
        ),
    ),
)
