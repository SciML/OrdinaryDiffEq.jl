import OrdinaryDiffEqSDC
using JET

@testset "JET Tests" begin
    test_package(OrdinaryDiffEqSDC, target_modules = (OrdinaryDiffEqSDC,), mode = :typo)
end
