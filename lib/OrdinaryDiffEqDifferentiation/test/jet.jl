using Pkg
Pkg.add("JET")

import OrdinaryDiffEqDifferentiation
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqDifferentiation, target_defined_modules = true, mode = :typo
    )
end
