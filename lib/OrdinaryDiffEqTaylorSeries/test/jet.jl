import OrdinaryDiffEqTaylorSeries
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqTaylorSeries, target_defined_modules = true, mode = :typo
    )
end
