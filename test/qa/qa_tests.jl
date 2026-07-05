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
    ),
)
