using SciMLTesting, OrdinaryDiffEq, SciMLBase
using Test

# The umbrella package's QA lane historically ran only the ExplicitImports checks
# (no Aqua/JET), so keep `aqua = false` to preserve that scope.
run_qa(
    OrdinaryDiffEq;
    aqua = false,
    explicit_imports = true,
    ei_kwargs = (;
        no_implicit_imports = (; skip = (Base, Core, SciMLBase)),
        # init/solve/solve!/step! (CommonSolve) and successful_retcode (SciMLBase)
        # are not declared public by their owners; OrdinaryDiffEq re-exports them.
        all_explicit_imports_are_public = (;
            ignore = (:init, :solve, :solve!, :step!, :successful_retcode),
        ),
    ),
)
