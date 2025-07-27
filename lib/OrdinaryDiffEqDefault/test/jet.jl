import OrdinaryDiffEqDefault
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqDefault, target_defined_modules = true, mode = :typo)
end
