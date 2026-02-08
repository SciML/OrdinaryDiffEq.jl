using Pkg
Pkg.add("JET")

import OrdinaryDiffEqNonlinearSolve
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqNonlinearSolve, mode = :typo
    )
end
