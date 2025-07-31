import OrdinaryDiffEqExplicitRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqExplicitRK, target_defined_modules = true, mode = :typo)
end
