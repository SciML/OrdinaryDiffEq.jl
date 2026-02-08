using Pkg
Pkg.add("JET")

import OrdinaryDiffEqFIRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqFIRK, target_modules = (OrdinaryDiffEqFIRK,), mode = :typo
    )
end
