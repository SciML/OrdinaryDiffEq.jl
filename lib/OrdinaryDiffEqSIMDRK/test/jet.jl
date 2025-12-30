import OrdinaryDiffEqSIMDRK
using OrdinaryDiffEqSIMDRK
using OrdinaryDiffEqCore
using JET
using Test

@testset "JET Tests" begin
    # Test package for typos
    test_package(
        OrdinaryDiffEqSIMDRK, target_modules = (OrdinaryDiffEqSIMDRK,), mode = :typo)
end
