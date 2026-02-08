using Pkg
Pkg.add("JET")

import OrdinaryDiffEqSymplecticRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqSymplecticRK, mode = :typo
    )
end
