using SciMLTesting, StochasticDiffEqIIF, Test
using JET

run_qa(
    StochasticDiffEqIIF;
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
    ei_kwargs = (;
        # These names are part of the OrdinaryDiffEqCore solver-author interface and
        # StochasticDiffEqCore's cache machinery, but are not (yet) declared `public`.
        # They are imported from their owning module; the only residual EI complaint is
        # that they are not public there. Tracked for make-public in
        # SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                :perform_step!, :issplit, :current_extrapolant, :current_extrapolant!,  # OrdinaryDiffEqCore
                Symbol("@cache"),  # StochasticDiffEqCore
            ),
        ),
    ),
)
