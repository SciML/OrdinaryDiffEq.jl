using Pkg
Pkg.add("JET")

import OrdinaryDiffEqRKN
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqRKN, mode = :typo
    )
end
