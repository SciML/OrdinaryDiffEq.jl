using Pkg
Pkg.add("JET")

import OrdinaryDiffEqCore
using JET, Test

@testset "JET Tests" begin
    @test test_package(
        OrdinaryDiffEqCore, target_defined_modules = true, mode = :typo
    ) === nothing broken = true
end
