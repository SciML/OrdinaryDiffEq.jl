using Pkg
Pkg.add("JET")

import OrdinaryDiffEqQPRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqQPRK, target_defined_modules = true, mode = :typo
    )
end
