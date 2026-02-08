using Pkg
Pkg.add("JET")

import OrdinaryDiffEqNonlinearSolve
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqNonlinearSolve, target_modules = (OrdinaryDiffEqNonlinearSolve,), mode = :typo
    )
end
