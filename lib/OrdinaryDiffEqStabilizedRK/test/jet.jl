import OrdinaryDiffEqStabilizedRK
using JET

@testset "JET Tests" begin
    test_package(test_package(
        OrdinaryDiffEqStabilizedRK, target_defined_modules = true, mode = :typo))
end
