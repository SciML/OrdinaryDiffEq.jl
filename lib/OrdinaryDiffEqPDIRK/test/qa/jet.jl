import OrdinaryDiffEqPDIRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqPDIRK, target_modules = (OrdinaryDiffEqPDIRK,), mode = :typo
    )
end
