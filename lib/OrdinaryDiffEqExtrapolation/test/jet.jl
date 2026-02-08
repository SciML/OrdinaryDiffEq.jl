using Pkg
Pkg.add("JET")

import OrdinaryDiffEqExtrapolation
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqExtrapolation, mode = :typo
    )
end
