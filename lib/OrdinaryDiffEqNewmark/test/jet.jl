import OrdinaryDiffEqNewmark
using OrdinaryDiffEqNewmark
using OrdinaryDiffEqCore
using JET
using Test

@testset "JET Tests" begin
    # Test package for typos
    test_package(
        OrdinaryDiffEqNewmark, target_modules = (OrdinaryDiffEqNewmark,), mode = :typo)
end
