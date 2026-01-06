import OrdinaryDiffEqHighOrderRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqHighOrderRK, target_defined_modules = true, mode = :typo
    )
end
