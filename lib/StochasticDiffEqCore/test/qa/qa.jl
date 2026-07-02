using StochasticDiffEqCore
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(StochasticDiffEqCore)
end

@testset "JET" begin
    JET.test_package(StochasticDiffEqCore; target_defined_modules = true)
end
