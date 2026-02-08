using Pkg
Pkg.add("JET")

import OrdinaryDiffEqTaylorSeries
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqTaylorSeries, mode = :typo
    )
end
