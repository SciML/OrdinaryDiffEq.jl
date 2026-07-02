using StochasticDiffEqLeaping
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(StochasticDiffEqLeaping)
end

@testset "JET" begin
    JET.test_package(StochasticDiffEqLeaping; target_defined_modules = true)
end
