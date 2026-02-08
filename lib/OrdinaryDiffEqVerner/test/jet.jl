using Pkg
Pkg.add("JET")

import OrdinaryDiffEqVerner
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqVerner, mode = :typo
    )
end
