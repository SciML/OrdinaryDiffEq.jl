import OrdinaryDiffEqAdamsBashforthMoulton
using JET

@testset "JET Tests" begin
    test_package(test_package(
        OrdinaryDiffEqAdamsBashforthMoulton, target_defined_modules = true, mode = :typo))
end