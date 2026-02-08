using Pkg
Pkg.add("JET")

import OrdinaryDiffEqLinear
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqLinear, mode = :typo
    )
end
