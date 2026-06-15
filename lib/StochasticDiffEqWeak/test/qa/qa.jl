using StochasticDiffEqWeak
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(StochasticDiffEqWeak)
end

@testset "JET" begin
    JET.test_package(StochasticDiffEqWeak; target_defined_modules = true)
end
