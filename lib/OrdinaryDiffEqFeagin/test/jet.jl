using Pkg
Pkg.add("JET")

import OrdinaryDiffEqFeagin
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqFeagin, mode = :typo
    )
end
