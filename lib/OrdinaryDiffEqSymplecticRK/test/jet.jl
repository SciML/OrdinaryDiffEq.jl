import OrdinaryDiffEqSymplecticRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqSymplecticRK, target_defined_modules = true, mode = :typo
    )
end
