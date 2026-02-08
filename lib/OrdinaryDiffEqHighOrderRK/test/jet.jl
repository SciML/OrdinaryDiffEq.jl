using Pkg
Pkg.add("JET")

import OrdinaryDiffEqHighOrderRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqHighOrderRK, target_modules = (OrdinaryDiffEqHighOrderRK,), mode = :typo
    )
end
