import OrdinaryDiffEqFIRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqFIRK, target_defined_modules = true, mode = :typo
    )
end
