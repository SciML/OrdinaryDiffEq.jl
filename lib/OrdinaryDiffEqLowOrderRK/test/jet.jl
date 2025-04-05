import OrdinaryDiffEqLowOrderRK
using JET

@testset "JET Tests" begin
    test_package(test_package(
        OrdinaryDiffEqLowOrderRK, target_defined_modules = true, mode = :typo))
end