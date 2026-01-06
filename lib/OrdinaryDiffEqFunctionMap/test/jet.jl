import OrdinaryDiffEqFunctionMap
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqFunctionMap, target_defined_modules = true, mode = :typo
    )
end
