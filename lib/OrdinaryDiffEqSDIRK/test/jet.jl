import OrdinaryDiffEqSDIRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqSDIRK, target_defined_modules = true, mode = :typo
    )
end
