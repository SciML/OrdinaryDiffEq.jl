import OrdinaryDiffEqRosenbrock
using JET

@testset "JET Tests" begin
    test_package(test_package(
        OrdinaryDiffEqRosenbrock, target_defined_modules = true, mode = :typo))
end
