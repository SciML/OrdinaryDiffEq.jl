using Pkg
Pkg.add("JET")

import OrdinaryDiffEqRosenbrock
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqRosenbrock, target_modules = (OrdinaryDiffEqRosenbrock,), mode = :typo
    )
end
