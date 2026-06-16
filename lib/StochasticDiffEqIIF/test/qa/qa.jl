using StochasticDiffEqIIF
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(StochasticDiffEqIIF)
end

@testset "JET" begin
    JET.test_package(StochasticDiffEqIIF; target_defined_modules = true)
end
