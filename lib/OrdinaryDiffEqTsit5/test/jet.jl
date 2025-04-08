import OrdinaryDiffEqTsit5
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqTsit5, target_defined_modules = true, mode = :typo)
end
