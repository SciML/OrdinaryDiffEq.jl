using Pkg
Pkg.add("JET")

import OrdinaryDiffEqStabilizedIRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqStabilizedIRK, target_modules = (OrdinaryDiffEqStabilizedIRK,), mode = :typo
    )
end
