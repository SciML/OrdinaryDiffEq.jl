using SciMLTesting, OrdinaryDiffEqLinear, Test

run_qa(
    OrdinaryDiffEqLinear;
    explicit_imports = true,
    ei_kwargs = (;
        # Residual non-public names below have no public-API replacement yet.
        # The OrdinaryDiffEqCore/DiffEqBase solver-author API is now declared
        # public on this branch, so those names are no longer ignored here.
        all_explicit_imports_are_public = (;
            ignore = (
                # SciMLOperators internal abstract type, non-public
                :AbstractSciMLOperator,
                # SciMLBase internal (owner of _vec), non-public
                :_vec,
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                # ExponentialUtilities internal cache allocator, non-public
                :alloc_mem,
            ),
        ),
    ),
)
