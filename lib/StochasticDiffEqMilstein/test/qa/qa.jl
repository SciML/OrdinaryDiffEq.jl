using StochasticDiffEqMilstein
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(StochasticDiffEqMilstein)
end

@testset "JET" begin
    JET.test_package(StochasticDiffEqMilstein; target_defined_modules = true)
end
