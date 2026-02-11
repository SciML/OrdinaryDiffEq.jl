using Pkg
Pkg.add("JET")

import OrdinaryDiffEqFeagin
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqFeagin, target_defined_modules = true, mode = :typo
    )
end
