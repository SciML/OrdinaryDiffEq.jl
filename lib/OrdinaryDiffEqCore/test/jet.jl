import OrdinaryDiffEqCore
using JET, Test

@testset "JET Tests" begin
    @test_broken test_package(
        OrdinaryDiffEqCore, target_defined_modules = true, mode = :typo) === nothing
end
