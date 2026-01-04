import OrdinaryDiffEqPRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqPRK, target_defined_modules = true, mode = :typo
    )
end
