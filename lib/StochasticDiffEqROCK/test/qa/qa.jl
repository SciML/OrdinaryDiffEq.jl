using StochasticDiffEqROCK
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(StochasticDiffEqROCK)
end

@testset "JET" begin
    JET.test_package(StochasticDiffEqROCK; target_defined_modules = true)
end
