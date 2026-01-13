using Pkg
Pkg.add("JET")

import OrdinaryDiffEqLowStorageRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqLowStorageRK, target_defined_modules = true, mode = :typo
    )
end
