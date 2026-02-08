using Pkg
Pkg.add("JET")

import OrdinaryDiffEqStabilizedRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqStabilizedRK, mode = :typo
    )
end
