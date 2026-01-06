import OrdinaryDiffEqVerner
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqVerner, target_defined_modules = true, mode = :typo
    )
end
