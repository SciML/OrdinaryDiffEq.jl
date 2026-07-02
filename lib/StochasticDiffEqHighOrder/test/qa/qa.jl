using SciMLTesting, StochasticDiffEqHighOrder, Test
using JET

run_qa(
    StochasticDiffEqHighOrder;
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
    ei_kwargs = (
        # `@..` is owned by FastBroadcast but reexported through DiffEqBase, and
        # FastBroadcast is not a direct dependency of this sublibrary, so the
        # broadcast macro can only be reached via the DiffEqBase reexport.
        all_explicit_imports_via_owners = (; ignore = (Symbol("@.."),)),
        # `@..`: external FastBroadcast macro reexported through DiffEqBase, where
        # it is not public.
        all_explicit_imports_are_public = (; ignore = (Symbol("@.."),)),
    ),
)
