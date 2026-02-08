using Pkg
Pkg.add("JET")

import OrdinaryDiffEqPDIRK
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqPDIRK, mode = :typo
    )
end
