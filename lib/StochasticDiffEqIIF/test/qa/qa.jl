using SciMLTesting, StochasticDiffEqIIF, Test
using JET

run_qa(
    StochasticDiffEqIIF;
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore solver-author interface, still not `public` there
                # (the no-bang `current_extrapolant` was made public; the in-place
                # `current_extrapolant!` was not).
                :current_extrapolant!,
                # `@cache` is public in OrdinaryDiffEqCore but imported here from
                # StochasticDiffEqCore, which does not re-declare it `public`.
                Symbol("@cache"),  # StochasticDiffEqCore
            ),
        ),
    ),
)
