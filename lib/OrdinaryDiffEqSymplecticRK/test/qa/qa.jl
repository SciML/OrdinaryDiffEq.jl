using SciMLTesting, OrdinaryDiffEqSymplecticRK, Test

run_qa(
    OrdinaryDiffEqSymplecticRK;
    explicit_imports = true,
)
