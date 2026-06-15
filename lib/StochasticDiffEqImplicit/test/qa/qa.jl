using StochasticDiffEqImplicit
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(StochasticDiffEqImplicit)
end

@testset "JET" begin
    JET.test_package(StochasticDiffEqImplicit; target_defined_modules = true)
end
