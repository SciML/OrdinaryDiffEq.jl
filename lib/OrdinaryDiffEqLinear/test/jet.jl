import OrdinaryDiffEqLinear
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqLinear, target_defined_modules = true, mode = :typo
    )
end
