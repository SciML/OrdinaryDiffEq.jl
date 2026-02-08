using Pkg
Pkg.add("JET")

import OrdinaryDiffEqQPRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqQPRK, target_modules = (OrdinaryDiffEqQPRK,), mode = :typo
    )
end
