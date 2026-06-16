using OrdinaryDiffEqRosenbrockTableaus
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(OrdinaryDiffEqRosenbrockTableaus)
end

@testset "JET" begin
    JET.test_package(OrdinaryDiffEqRosenbrockTableaus; target_defined_modules = true)
end
