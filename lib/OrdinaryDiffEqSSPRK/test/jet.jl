import OrdinaryDiffEqSSPRK
using JET

@testset "JET Tests" begin
    test_package(test_package(
        OrdinaryDiffEqSSPRK, target_defined_modules = true, mode = :typo))
end
