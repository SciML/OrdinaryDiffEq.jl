import OrdinaryDiffEqPDIRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqPDIRK, target_defined_modules = true, mode = :typo
    )
end
