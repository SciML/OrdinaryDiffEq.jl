using Pkg
Pkg.add("JET")

import OrdinaryDiffEqLinear
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqLinear, target_modules = (OrdinaryDiffEqLinear,), mode = :typo
    )
end
