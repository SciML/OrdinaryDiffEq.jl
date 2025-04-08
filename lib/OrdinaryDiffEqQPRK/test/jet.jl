import OrdinaryDiffEqQRPK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqQRPK, target_defined_modules = true, mode = :typo)
end
