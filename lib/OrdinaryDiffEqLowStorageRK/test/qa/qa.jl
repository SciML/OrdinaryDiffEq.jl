using SciMLTesting, OrdinaryDiffEqLowStorageRK, Test

run_qa(
    OrdinaryDiffEqLowStorageRK;
    explicit_imports = true,
    ei_kwargs = (
        # Every remaining name is a genuine non-public internal of its owner. The
        # solver-author API of OrdinaryDiffEqCore / DiffEqBase is now declared
        # `public`, so those ignores were dropped.
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore precompile-workload probes (owner-internal)
                :lorenz, :lorenz_oop,
                # Base.Broadcast internals used in ArrayFuse copyto!/materialize!
                :Broadcasted, :materialize!,
            )),
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore default no-op limiter (owner-internal)
                :trivial_limiter!,
            )),
    ),
)
