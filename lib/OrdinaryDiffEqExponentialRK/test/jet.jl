import OrdinaryDiffEqExponentialRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqExponentialRK, target_defined_modules = true, mode = :typo
    )
end
