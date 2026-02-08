using Pkg
Pkg.add("JET")

import OrdinaryDiffEqFunctionMap
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqFunctionMap, mode = :typo
    )
end
