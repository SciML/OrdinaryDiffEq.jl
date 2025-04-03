import OrdinaryDiffEqExponentialRK
using JET

@testset "JET Tests" begin
    test_package(test_package(
        OrdinaryDiffEqDifferentiation, target_defined_modules = true, mode = :typo))
end