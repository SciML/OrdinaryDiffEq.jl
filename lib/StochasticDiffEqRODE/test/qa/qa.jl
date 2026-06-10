using StochasticDiffEqRODE
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(StochasticDiffEqRODE)
end

@testset "JET" begin
    JET.test_package(StochasticDiffEqRODE; target_defined_modules = true)
end
