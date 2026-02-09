using Pkg
Pkg.add("JET")

import OrdinaryDiffEqSymplecticRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqSymplecticRK, target_modules = (OrdinaryDiffEqSymplecticRK,), mode = :typo
    )
end
