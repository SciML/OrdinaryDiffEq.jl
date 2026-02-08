using Pkg
Pkg.add("JET")

import OrdinaryDiffEqNordsieck
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqNordsieck, mode = :typo
    )
end
