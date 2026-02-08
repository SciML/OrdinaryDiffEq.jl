using Pkg
Pkg.add("JET")

import OrdinaryDiffEqSDIRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqSDIRK, target_modules = (OrdinaryDiffEqSDIRK,), mode = :typo
    )
end
