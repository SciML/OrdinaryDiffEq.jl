using Pkg
Pkg.add("JET")

import OrdinaryDiffEqExponentialRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqExponentialRK, mode = :typo
    )
end
