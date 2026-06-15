using StochasticDiffEqLevyArea
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(StochasticDiffEqLevyArea)
end

@testset "JET" begin
    JET.test_package(StochasticDiffEqLevyArea; target_defined_modules = true)
end
