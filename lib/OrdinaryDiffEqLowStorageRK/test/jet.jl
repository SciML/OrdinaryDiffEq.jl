using Pkg
Pkg.add("JET")

import OrdinaryDiffEqLowStorageRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqLowStorageRK, target_modules = (OrdinaryDiffEqLowStorageRK,), mode = :typo
    )
end
