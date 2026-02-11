using Pkg
Pkg.add("JET")

import OrdinaryDiffEqNordsieck
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqNordsieck, target_defined_modules = true, mode = :typo
    )
end
