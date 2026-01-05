import OrdinaryDiffEqNonlinearSolve
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqNonlinearSolve, target_defined_modules = true, mode = :typo
    )
end
