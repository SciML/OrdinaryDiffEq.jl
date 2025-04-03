import OrdinaryDiffEqCore
using JET

@testset "JET Tests" begin
    test_package(test_package(
        OrdinaryDiffEqCore, target_defined_modules = true, mode = :typo))
end