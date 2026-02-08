using Pkg
Pkg.add("JET")

import OrdinaryDiffEqQPRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqQPRK, mode = :typo
    )
end
