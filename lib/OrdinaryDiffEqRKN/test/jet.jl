import OrdinaryDiffEqRKN
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqRKN, target_defined_modules = true, mode = :typo
    )
end
