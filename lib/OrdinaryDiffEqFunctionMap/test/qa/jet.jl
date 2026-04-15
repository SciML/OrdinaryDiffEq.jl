import OrdinaryDiffEqFunctionMap
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqFunctionMap, target_modules = (OrdinaryDiffEqFunctionMap,), mode = :typo
    )
end
