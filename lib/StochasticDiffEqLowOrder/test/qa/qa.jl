using StochasticDiffEqLowOrder
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(StochasticDiffEqLowOrder)
end

@testset "JET" begin
    JET.test_package(StochasticDiffEqLowOrder; target_defined_modules = true)
end
