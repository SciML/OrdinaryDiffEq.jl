import OrdinaryDiffEqRKIP
using OrdinaryDiffEqRKIP
using OrdinaryDiffEqCore
using JET
using Test

@testset "JET Tests" begin
    # Test package for typos
    test_package(
        OrdinaryDiffEqRKIP, target_modules = (OrdinaryDiffEqRKIP,), mode = :typo
    )
end
