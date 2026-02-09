using Pkg
Pkg.add("JET")

import OrdinaryDiffEqExtrapolation
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqExtrapolation, target_modules = (OrdinaryDiffEqExtrapolation,), mode = :typo
    )
end
