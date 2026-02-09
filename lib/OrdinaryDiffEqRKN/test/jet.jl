using Pkg
Pkg.add("JET")

import OrdinaryDiffEqRKN
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqRKN, target_modules = (OrdinaryDiffEqRKN,), mode = :typo
    )
end
