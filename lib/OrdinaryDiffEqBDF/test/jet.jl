import OrdinaryDiffEqBDF
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqBDF, target_defined_modules = true, mode = :typo)
end
