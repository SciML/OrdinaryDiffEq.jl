using StochasticDiffEqHighOrder
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(StochasticDiffEqHighOrder)
end

@testset "JET" begin
    JET.test_package(StochasticDiffEqHighOrder; target_defined_modules = true)
end
