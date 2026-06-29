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
        # Solver-internal names with no public counterpart in their owner.
        # DiffEqBase residual hooks (calculate_residuals[!]), the Tableau type,
        # and the @.. broadcast macro; StochasticDiffEqCore's @cache macro; and
        # OrdinaryDiffEqCore's perform_step! dispatch hook.
        # See SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                Symbol("@.."), Symbol("@cache"),
                :Tableau, :calculate_residuals, :calculate_residuals!,
                :perform_step!,
            ),
        ),
        # Qualified access to OrdinaryDiffEqCore's internal step-error setter.
        all_qualified_accesses_are_public = (; ignore = (:set_EEst!,)),
    ),
)
