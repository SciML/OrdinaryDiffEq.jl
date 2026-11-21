using SciMLTesting, OrdinaryDiffEqRosenbrockTableaus, Test
using JET

run_qa(
    OrdinaryDiffEqRosenbrockTableaus;
    explicit_imports = true,
)
